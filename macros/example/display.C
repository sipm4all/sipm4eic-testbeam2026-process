#include "../lib/frame_selection.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TPad.h>
#include <TStyle.h>
#include <TText.h>
#include <TEllipse.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

constexpr double time_min = -100.;
constexpr double time_max = 100.;
constexpr double pixel_size = 3.0;
// These must match the ring-finder-ransac settings used to produce the ring tree.
constexpr double ring_match_tolerance = 5.;
constexpr double ring_match_time_window = 4.;
constexpr bool draw_matched_ring_hits = true;


namespace {

bool
finite_xy(const hit_t &hit)
{
  return std::isfinite(hit.x) && std::isfinite(hit.y);
}

double
ring_residual(const ring_t &ring, const hit_t &hit)
{
  const double dx = hit.x - ring.x0;
  const double dy = hit.y - ring.y0;
  const double c = std::cos(ring.phi);
  const double s = std::sin(ring.phi);
  const double xp = c * dx + s * dy;
  const double yp = -s * dx + c * dy;
  const double minor = ring.radius * std::sqrt(std::max(
    1.e-12, 1. - ring.eccentricity * ring.eccentricity));
  if (ring.radius <= 0. || minor <= 0.)
    return std::numeric_limits<double>::infinity();
  return std::abs(std::hypot(xp / ring.radius, yp / minor) - 1.) * ring.radius;
}

bool
is_ring_hit(const hit_t &hit, const ring_t &ring)
{
  return finite_xy(hit) && std::isfinite(hit.time) &&
         ring_residual(ring, hit) <= ring_match_tolerance &&
         std::abs(hit.time - ring.time) <= ring_match_time_window;
}

bool
passes_trigger_requirement(const trigger_reader_t &reader,
                           const trigger_selection_t &selection)
{
  if (!selection.enabled)
    return true;

  for (const auto &hit : reader.trigger_hits()) {
    if (selection.device < 0 || hit.device == selection.device)
      return true;
  }
  return false;
}

void
save_frame_png(TCanvas *canvas, const trigger_reader_t &reader)
{
  const std::string filename = Form("spill_%d_frame_%d.png",
                                    reader.spill_id(), reader.frame_index());
  canvas->SaveAs(filename.c_str());
  std::cout << "saved " << filename << std::endl;
}

bool
passes_ring_selection(const trigger_reader_t &reader,
                       const ring_selection_t &selection)
{
  const int n = static_cast<int>(reader.rings().size());
  if (n < selection.min_n_rings || n > selection.max_n_rings)
    return false;
  if (n == 0)
    return true;

  for (const auto &ring : reader.rings()) {
    if (ring.x0 >= selection.min_x0 && ring.x0 <= selection.max_x0 &&
        ring.y0 >= selection.min_y0 && ring.y0 <= selection.max_y0 &&
        ring.radius >= selection.min_r && ring.radius <= selection.max_r)
      return true;
  }
  return false;
}

int
color_index(double value, double vmin, double vmax)
{
  constexpr int ncolors = 255;
  if (vmax <= vmin)
    return gStyle->GetColorPalette(ncolors / 2);

  auto f = (value - vmin) / (vmax - vmin);
  f = std::max(0., std::min(1., f));
  auto index = static_cast<int>(std::lround(f * (ncolors - 1)));
  return gStyle->GetColorPalette(index);
}

bool
match_field(int value, int requested)
{
  return requested < 0 || value == requested;
}

int
channel(const hit_t &hit)
{
  int chip = hit.fifo / 4;
  int local_channel = hit.pixel + 4 * hit.column + 32 * chip;
  return local_channel + 256 * (hit.device - 192);
}

bool
match(const hit_t &hit, field_selector_t selector)
{
  return match_field(hit.type, selector.type) &&
         match_field(hit.device, selector.device) &&
         match_field(hit.fifo, selector.fifo) &&
         match_field(hit.column, selector.column) &&
         match_field(hit.pixel, selector.pixel);
}

bool
match(const hit_t &hit, channel_selector_t selector)
{
  if (selector.type == 9)
    return match(hit, field_selector_t{9, selector.channel, 32, -1, -1});

  if (!match_field(hit.type, selector.type))
    return false;

  if (hit.type != 1)
    return false;

  return match_field(channel(hit), selector.channel);
}

template <typename selector_t>
bool
find_reference_time(const trigger_reader_t &reader,
                    selector_t selector,
                    double &reference_time)
{
  auto scan = [&](const std::vector<hit_t> &hits) {
    for (const auto &hit : hits) {
      if (std::isfinite(hit.time) && match(hit, selector)) {
        reference_time = hit.time;
        return true;
      }
    }
    return false;
  };

  return scan(reader.trigger_hits()) ||
         scan(reader.timing_hits()) ||
         scan(reader.cherenkov_hits());
}

bool
find_reference_time(const trigger_reader_t &reader,
                    timing_reference_t selection,
                    double &reference_time)
{
  if (!reader.has_timing()) {
    std::cerr << "ERROR: input file does not contain timing estimator branches" << std::endl;
    return false;
  }

  if (reader.timing_valid() == 0)
    return false;

  if (selection.branch == "T0")
    reference_time = reader.T0();
  else if (selection.branch == "T1")
    reference_time = reader.T1();
  else if (selection.branch == "T")
    reference_time = reader.T();
  else {
    std::cerr << "ERROR: unsupported timing selection branch: " << selection.branch
              << " (valid choices are T, T0, T1)" << std::endl;
    return false;
  }

  return std::isfinite(reference_time);
}

void
draw_frame_delta(const trigger_reader_t &reader,
                 bool use_reference,
                 double reference_time)
{
  constexpr double time_to_ns = 3.125;
  static TCanvas *canvas = nullptr;
  if (!canvas) {
    canvas = new TCanvas("cDelta", "Cherenkov #Deltat", 800, 800);
    canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
    canvas->SetFillColor(kWhite);
    canvas->SetFrameFillColor(kWhite);
  }
  canvas->cd();

  auto histogram = new TH1D("display_delta", "", 400, -100, 100);
  histogram->SetDirectory(nullptr);
  histogram->SetStats(0);
  histogram->SetTitle(Form("spill %d frame %d Cherenkov timing;%s [ns];hits",
                           reader.spill_id(), reader.frame_index(),
                           use_reference ? "#Deltat" : "time"));
  histogram->SetLineColor(kBlack);
  histogram->SetFillColor(kAzure - 9);

  for (const auto &hit : reader.cherenkov_hits()) {
    if (!std::isfinite(hit.time))
      continue;
    const double value = use_reference ? hit.time - reference_time : hit.time;
    histogram->Fill(value * time_to_ns);
  }

  gPad->Clear();
  histogram->Draw();
  gPad->Modified();
  gPad->Update();
}

void
draw_frame_map(trigger_reader_t &reader,
               bool use_reference,
               double reference_time)
{
  constexpr double time_to_ns = 3.125;
  auto hits = reader.cherenkov_hits();
  std::vector<hit_t> drawable;
  drawable.reserve(hits.size());

  double vmin = std::numeric_limits<double>::infinity();
  double vmax = -std::numeric_limits<double>::infinity();

  for (const auto &hit : hits) {
    if (!finite_xy(hit) || !std::isfinite(hit.time))
      continue;
    drawable.push_back(hit);
    vmin = std::min(vmin, hit.time);
    vmax = std::max(vmax, hit.time);
  }

  constexpr double color_min = time_min;
  constexpr double color_max = time_max;

  gStyle->SetPalette(kBird);
  gStyle->SetPalette(kRainbow);

  gPad->Clear();

  if (drawable.empty()) {
    auto frame = gPad->DrawFrame(-99., -99., 99., 99.);
    frame->SetTitle(Form("spill %d frame %d; x (mm); y (mm)",
                         reader.spill_id(), reader.frame_index()));
    frame->GetXaxis()->SetTitleOffset(1.5);
    frame->GetYaxis()->SetTitleOffset(1.5);
    gPad->Update();
    draw_frame_delta(reader, use_reference, reference_time);
    return;
  }

  auto frame = gPad->DrawFrame(-99., -99., 99., 99.);
  frame->GetXaxis()->SetTitleOffset(1.5);
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetTitle(Form("spill %d frame %d; x (mm); y (mm)",
                       reader.spill_id(), reader.frame_index()));
  gPad->SetBit(TPad::kClipFrame);

  // PDU placement coordinates are the bottom-left corners, in millimetres.
  constexpr double pdu_size = 52.;
  const std::array<std::pair<int, std::array<double, 2>>, 8> pdus = {{
    {1, {-82.,  30.}}, {2, {-26.,  35.}}, {3, {30.,  30.}},
    {8, {-82., -26.}},                    {4, {30., -26.}},
    {7, {-82., -82.}}, {6, {-26., -87.}}, {5, {30., -82.}}
  }};
  for (const auto &pdu : pdus) {
    const double x = pdu.second[0];
    const double y = pdu.second[1];
    auto box = new TBox(x, y, x + pdu_size, y + pdu_size);
    box->SetFillStyle(0);
    box->SetLineColor(kGray + 1);
    box->SetLineWidth(1);
    box->Draw("same");
  }

  // Draw the colour scale directly, without using a filled helper histogram.
  constexpr int n_palette_bins = 48;
  const double palette_x1 = 105;
  const double palette_x2 = 105. + 10.;
  const double palette_y1 = -99.;
  const double palette_y2 = 99.;
  for (int i = 0; i < n_palette_bins; ++i) {
    const double y1 = palette_y1 + (palette_y2 - palette_y1) * i / n_palette_bins;
    const double y2 = palette_y1 + (palette_y2 - palette_y1) * (i + 1) / n_palette_bins;
    const double value = color_min + (color_max - color_min) * (i + 0.5) / n_palette_bins;
    auto box = new TBox(palette_x1, y1, palette_x2, y2);
    box->SetFillColor(color_index(value, color_min, color_max));
    box->SetLineColor(color_index(value, color_min, color_max));
    box->Draw("same");
  }
  auto palette_title = new TText((palette_x1 + palette_x2) * 0.5, palette_y2 + 15.,
                                 use_reference ? "t (ns)" : "time (ns)");
  palette_title->SetTextSize(0.035);
  palette_title->SetTextFont(42);
  palette_title->SetTextAlign(21);
  palette_title->Draw("same");
  for (int i = 0; i <= 4; ++i) {
    const double f = i / 4.;
    auto label = new TText(palette_x2 + 2.,
                           palette_y1 + f * (palette_y2 - palette_y1),
                           Form("%.3g", color_min + f * (color_max - color_min)));
    label->SetTextSize(0.035);
    label->SetTextFont(42);
    label->SetTextAlign(12);
    label->Draw("same");
  }

  for (const auto &hit : drawable) {
    auto dt = (use_reference ? hit.time - reference_time : hit.time) * time_to_ns;
    const int hit_color = color_index(dt, color_min, color_max);
    auto box = new TBox(hit.x - 0.5 * pixel_size, hit.y - 0.5 * pixel_size,
                        hit.x + 0.5 * pixel_size, hit.y + 0.5 * pixel_size);
    box->SetFillColor(hit_color);
    bool ring_hit = false;
    for (const auto &ring : reader.rings()) {
      if (is_ring_hit(hit, ring)) {
        ring_hit = true;
        break;
      }
    }
    box->SetLineColor(kBlack);
    box->SetLineWidth(1);
    box->Draw("same");

    if (draw_matched_ring_hits && ring_hit) {
      // Empty circle is slightly larger than the square and retains its time color.
      auto marker = new TEllipse(hit.x, hit.y, 1.0 * pixel_size, 1.0 * pixel_size);
      marker->SetFillStyle(0);
      marker->SetLineColor(hit_color);
      marker->SetLineWidth(2);
      marker->Draw("same");
    }
  }

  if (reader.has_rings()) {
    for (const auto &ring : reader.rings()) {
      double minor = ring.radius * std::sqrt(std::max(
        0., 1. - ring.eccentricity * ring.eccentricity));
      auto circle = new TEllipse(ring.x0, ring.y0, ring.radius, minor,
                                 0., 360., ring.phi * 180. / std::acos(-1.));
      circle->SetFillStyle(0);
      const double ring_value = use_reference ?
        (ring.time - reference_time) * time_to_ns :
        ring.time * time_to_ns;
      circle->SetLineColor(color_index(ring_value, color_min, color_max));
      circle->SetLineWidth(2);
      circle->Draw("same");

    }
  }

  auto label = new TText(-98., 94.,
                         use_reference ?
                         Form("hits=%zu  dt=[%.3f, %.3f]",
                              drawable.size(),
                              (vmin - reference_time) * time_to_ns,
                              (vmax - reference_time) * time_to_ns) :
                         Form("hits=%zu  time range=%.3f  rings=%zu",
                              drawable.size(), (vmax - vmin) * time_to_ns,
                              reader.rings().size()));
  label->SetTextSize(0.03);
  //  label->Draw("same");

  gPad->Modified();
  gPad->Update();
  draw_frame_delta(reader, use_reference, reference_time);
}

} // namespace


namespace {

TCanvas *
make_display_canvas()
{
  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFrameFillColor(kWhite);
  canvas->cd();
  gPad->SetFillColor(kWhite);
  return canvas;
}

bool
passes_selections(const trigger_reader_t &reader,
                  const std::vector<selection_ptr_t> &selections)
{
  for (const auto &selection : selections) {
    if (!selection) {
      std::cerr << "ERROR: null display selection" << std::endl;
      return false;
    }
    if (!selection->is_selected(reader))
      return false;
  }
  return true;
}

bool
validate_selections(const trigger_reader_t &reader,
                    const std::vector<selection_ptr_t> &selections)
{
  for (const auto &selection : selections) {
    if (const auto trigger = dynamic_cast<const trigger_selection_t *>(selection.get());
        trigger && trigger->enabled && !reader.has_trigger_hits()) {
      std::cerr << "ERROR: trigger selection requested, but input has no trigger tree"
                << std::endl;
      return false;
    }
  }
  return true;
}

void
display_frames(const char *filename,
               const reference_ptr_t &reference,
               const std::vector<selection_ptr_t> &selections,
               int start_spill,
               int start_frame,
               int target_spill,
               int target_frame,
               const char *ring_name)
{
  trigger_reader_t reader;
  if (!reader.open(filename, ring_name ? ring_name : "ring"))
    return;
  if (!validate_selections(reader, selections))
    return;

  auto canvas = make_display_canvas();
  const bool fixed_frame = target_spill != std::numeric_limits<int>::min();

  while (reader.next_spill()) {
    if (reader.spill_id() < start_spill)
      continue;
    if (fixed_frame && reader.spill_id() != target_spill)
      continue;

    while (reader.next_frame()) {
      if (reader.spill_id() == start_spill &&
          reader.frame_index() < start_frame)
        continue;
      if (fixed_frame && reader.frame_index() != target_frame)
        continue;
      if (!passes_selections(reader, selections))
        continue;

      bool use_reference = false;
      double reference_time = 0.;
      if (reference) {
        if (!reference->process(reader)) {
          std::cerr << "spill " << reader.spill_id()
                    << " frame " << reader.frame_index()
                    << " -- no matching reference, skipping"
                    << std::endl;
          continue;
        }
        use_reference = true;
        reference_time = reference->reference_time();
      }

      canvas->cd();
      draw_frame_map(reader, use_reference, reference_time);

      if (fixed_frame)
        return;

      std::cout << "spill " << reader.spill_id()
                << " frame " << reader.frame_index();
      if (use_reference)
        std::cout << " reference_time=" << reference_time;
      std::cout << " trigger=" << reader.trigger_hits().size()
                << " timing=" << reader.timing_hits().size()
                << " cherenkov=" << reader.cherenkov_hits().size()
                << " -- Enter: next, s: save PNG, q: quit"
                << std::endl;

      std::string line;
      std::getline(std::cin, line);
      if (line == "q" || line == "Q")
        return;
      if (line == "s" || line == "S")
        save_frame_png(canvas, reader);
    }
  }

  if (fixed_frame) {
    std::cerr << "ERROR: frame not found"
              << " spill=" << target_spill
              << " frame=" << target_frame
              << std::endl;
  }
}

} // namespace

void
display(const char *filename = "triggered.root",
        reference_ptr_t reference = nullptr,
        const std::vector<selection_ptr_t> &selections = {},
        int start_spill = std::numeric_limits<int>::min(),
        int start_frame = 0,
        int target_spill = std::numeric_limits<int>::min(),
        int target_frame = -1,
        const char *ring_name = "ring")
{
  display_frames(filename, reference, selections, start_spill, start_frame,
                 target_spill, target_frame, ring_name);
}

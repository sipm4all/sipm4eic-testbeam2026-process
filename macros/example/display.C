#include "../lib/trigger_reader.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TPad.h>
#include <TStyle.h>
#include <TText.h>
#include <TEllipse.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct field_selector_t {
  int type;
  int device;
  int fifo;
  int column;
  int pixel;

  field_selector_t(int type_, int device_, int fifo_, int column_, int pixel_)
    : type(type_), device(device_), fifo(fifo_), column(column_), pixel(pixel_)
  {
  }
};

struct channel_selector_t {
  int type;
  int channel;

  channel_selector_t(int type_, int channel_)
    : type(type_), channel(channel_)
  {
  }
};

struct timing_selection_t {
  std::string branch;

  timing_selection_t(const std::string &branch_ = "T")
    : branch(branch_)
  {
  }
};

struct reference_time_t {
  std::string branch;

  reference_time_t(const std::string &branch_ = "T")
    : branch(branch_)
  {
  }
};

struct ring_selection_t {
  double min_x0 = -std::numeric_limits<double>::infinity();
  double max_x0 =  std::numeric_limits<double>::infinity();
  double min_y0 = -std::numeric_limits<double>::infinity();
  double max_y0 =  std::numeric_limits<double>::infinity();
  double min_r = 0.;
  double max_r = std::numeric_limits<double>::infinity();
  int min_n_rings = 1;
  int max_n_rings = std::numeric_limits<int>::max();

  ring_selection_t(double min_x0_ = -std::numeric_limits<double>::infinity(),
                   double max_x0_ =  std::numeric_limits<double>::infinity(),
                   double min_y0_ = -std::numeric_limits<double>::infinity(),
                   double max_y0_ =  std::numeric_limits<double>::infinity(),
                   double min_r_ = 0.,
                   double max_r_ = std::numeric_limits<double>::infinity(),
                   int min_n_rings_ = 1,
                   int max_n_rings_ = std::numeric_limits<int>::max())
    : min_x0(min_x0_), max_x0(max_x0_), min_y0(min_y0_), max_y0(max_y0_),
      min_r(min_r_), max_r(max_r_), min_n_rings(min_n_rings_),
      max_n_rings(max_n_rings_)
  {
  }
};

namespace {

bool
finite_xy(const hit_t &hit)
{
  return std::isfinite(hit.x) && std::isfinite(hit.y);
}

bool
passes_ring_selection(const trigger_reader_t &reader,
                       const ring_selection_t &selection)
{
  const int n = static_cast<int>(reader.rings().size());
  if (n < selection.min_n_rings || n > selection.max_n_rings)
    return false;

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
                    timing_selection_t selection,
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

bool
find_reference_time(const trigger_reader_t &reader,
                    reference_time_t selection,
                    double &reference_time)
{
  return find_reference_time(reader, timing_selection_t(selection.branch),
                             reference_time);
}

void
draw_frame_map(trigger_reader_t &reader,
               bool use_reference,
               double reference_time,
               double pixel_size)
{
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

  constexpr double color_min = -15.;
  constexpr double color_max = 15.;

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
                                 use_reference ? "t (au)" : "time");
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
    auto dt = use_reference ? hit.time - reference_time : hit.time - vmin;
    auto box = new TBox(hit.x - 0.5 * pixel_size, hit.y - 0.5 * pixel_size,
                        hit.x + 0.5 * pixel_size, hit.y + 0.5 * pixel_size);
    box->SetFillColor(color_index(dt, color_min, color_max));
    box->SetLineColor(kBlack);
    box->SetLineWidth(1);
    box->Draw("same");
  }

  if (reader.has_rings()) {
    for (const auto &ring : reader.rings()) {
      double minor = ring.radius * std::sqrt(std::max(
        0., 1. - ring.eccentricity * ring.eccentricity));
      auto circle = new TEllipse(ring.x0, ring.y0, ring.radius, minor,
                                 0., 360., ring.phi * 180. / std::acos(-1.));
      circle->SetFillStyle(0);
      const double ring_value = use_reference ?
        ring.time - reference_time : ring.time - vmin;
      circle->SetLineColor(color_index(ring_value, color_min, color_max));
      circle->SetLineWidth(2);
      circle->Draw("same");

    }
  }

  auto label = new TText(-98., 94.,
                         use_reference ?
                         Form("hits=%zu  dt=[%.3f, %.3f]",
                              drawable.size(), vmin, vmax) :
                         Form("hits=%zu  time range=%.3f  rings=%zu",
                              drawable.size(), vmax - vmin, reader.rings().size()));
  label->SetTextSize(0.03);
  //  label->Draw("same");

  gPad->Modified();
  gPad->Update();
}

} // namespace

void
display(const char *filename = "triggered.root",
        double pixel_size = 3.2)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFrameFillColor(kWhite);
  canvas->cd();
  gPad->SetFillColor(kWhite);

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      canvas->cd();
      draw_frame_map(reader, false, 0., pixel_size);

      std::cout << "spill " << reader.spill_id()
                << " frame " << reader.frame_index()
                << " trigger=" << reader.trigger_hits().size()
                << " timing=" << reader.timing_hits().size()
                << " cherenkov=" << reader.cherenkov_hits().size()
                << " -- press Enter for next frame, q then Enter to quit"
                << std::endl;

      std::string line;
      std::getline(std::cin, line);
      if (line == "q" || line == "Q")
        return;
    }
  }
}

void
display(const char *filename,
        int target_spill,
        int target_frame,
        double pixel_size = 3.2)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFrameFillColor(kWhite);
  canvas->cd();
  gPad->SetFillColor(kWhite);

  while (reader.next_spill()) {
    if (reader.spill_id() != target_spill)
      continue;

    while (reader.next_frame()) {
      if (reader.frame_index() != target_frame)
        continue;

      canvas->cd();
      draw_frame_map(reader, false, 0., pixel_size);
      return;
    }
  }

  std::cerr << "ERROR: frame not found"
            << " spill=" << target_spill
            << " frame=" << target_frame
            << std::endl;
}

template <typename selector_t>
void
display_reference_loop(const char *filename,
                       selector_t reference,
                       double pixel_size,
                       int start_spill = std::numeric_limits<int>::min(),
                       int start_frame = 0,
                       const ring_selection_t *ring_selection = nullptr)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFrameFillColor(kWhite);
  canvas->cd();
  gPad->SetFillColor(kWhite);

  while (reader.next_spill()) {
    if (reader.spill_id() < start_spill)
      continue;
    while (reader.next_frame()) {
      if (reader.spill_id() == start_spill &&
          reader.frame_index() < start_frame)
        continue;
      if (ring_selection &&
          (!reader.has_rings() ||
           !passes_ring_selection(reader, *ring_selection)))
        continue;
      double reference_time = 0.;
      if (!find_reference_time(reader, reference, reference_time)) {
        std::cerr << "spill " << reader.spill_id()
                  << " frame " << reader.frame_index()
                  << " -- no matching reference hit, skipping"
                  << std::endl;
        continue;
      }

      canvas->cd();
      draw_frame_map(reader, true, reference_time, pixel_size);

      std::cout << "spill " << reader.spill_id()
                << " frame " << reader.frame_index()
                << " reference_time=" << reference_time
                << " trigger=" << reader.trigger_hits().size()
                << " timing=" << reader.timing_hits().size()
                << " cherenkov=" << reader.cherenkov_hits().size()
                << " -- press Enter for next frame, q then Enter to quit"
                << std::endl;

      std::string line;
      std::getline(std::cin, line);
      if (line == "q" || line == "Q")
        return;
    }
  }
}

template <typename selector_t>
void
display_reference_frame(const char *filename,
                        int target_spill,
                        int target_frame,
                        selector_t reference,
                        double pixel_size)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFrameFillColor(kWhite);
  canvas->cd();
  gPad->SetFillColor(kWhite);

  while (reader.next_spill()) {
    if (reader.spill_id() != target_spill)
      continue;

    while (reader.next_frame()) {
      if (reader.frame_index() != target_frame)
        continue;

      double reference_time = 0.;
      if (!find_reference_time(reader, reference, reference_time)) {
        std::cerr << "ERROR: no matching reference hit"
                  << " spill=" << target_spill
                  << " frame=" << target_frame
                  << std::endl;
        return;
      }

      canvas->cd();
      draw_frame_map(reader, true, reference_time, pixel_size);
      return;
    }
  }

  std::cerr << "ERROR: frame not found"
            << " spill=" << target_spill
            << " frame=" << target_frame
            << std::endl;
}

void
display(const char *filename,
        field_selector_t reference,
        double pixel_size = 3.2)
{
  display_reference_loop(filename, reference, pixel_size);
}

void
display(const char *filename,
        channel_selector_t reference,
        double pixel_size = 3.2)
{
  display_reference_loop(filename, reference, pixel_size);
}

void
display(const char *filename,
        timing_selection_t reference,
        double pixel_size = 3.2)
{
  display_reference_loop(filename, reference, pixel_size);
}

void
display(const char *filename,
        reference_time_t reference,
        int start_spill,
        int start_frame,
        double pixel_size = 3.2)
{
  display_reference_loop(filename, reference, pixel_size,
                         start_spill, start_frame);
}

void
display(const char *filename,
        reference_time_t reference,
        ring_selection_t ring_selection,
        double pixel_size = 3.2)
{
  display_reference_loop(filename, reference, pixel_size,
                         std::numeric_limits<int>::min(), 0,
                         &ring_selection);
}

void
display(const char *filename,
        ring_selection_t ring_selection,
        double pixel_size = 3.2)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFrameFillColor(kWhite);
  canvas->cd();
  gPad->SetFillColor(kWhite);

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      if (!reader.has_rings() ||
          !passes_ring_selection(reader, ring_selection))
        continue;
      canvas->cd();
      draw_frame_map(reader, false, 0., pixel_size);
      std::cout << "spill " << reader.spill_id()
                << " frame " << reader.frame_index()
                << " -- press Enter for next frame, q then Enter to quit"
                << std::endl;
      std::string line;
      std::getline(std::cin, line);
      if (line == "q" || line == "Q")
        return;
    }
  }
}

void
display(const char *filename,
        int target_spill,
        int target_frame,
        field_selector_t reference,
        double pixel_size = 3.2)
{
  display_reference_frame(filename, target_spill, target_frame, reference, pixel_size);
}

void
display(const char *filename,
        int target_spill,
        int target_frame,
        timing_selection_t reference,
        double pixel_size = 3.2)
{
  display_reference_frame(filename, target_spill, target_frame, reference, pixel_size);
}

void
display(const char *filename,
        int target_spill,
        int target_frame,
        channel_selector_t reference,
        double pixel_size = 3.2)
{
  display_reference_frame(filename, target_spill, target_frame, reference, pixel_size);
}

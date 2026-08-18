#include "../lib/trigger_reader.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TText.h>

#include <algorithm>
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

namespace {

bool
finite_xy(const hit_t &hit)
{
  return std::isfinite(hit.x) && std::isfinite(hit.y);
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
    double value = use_reference ? hit.time - reference_time : hit.time;
    vmin = std::min(vmin, value);
    vmax = std::max(vmax, value);
  }

  double color_min = use_reference ? vmin : 0.;
  double color_max = use_reference ? vmax : vmax - vmin;

  gStyle->SetPalette(kViridis);

  gPad->Clear();
  gPad->SetRightMargin(0.16);

  if (drawable.empty()) {
    auto empty = new TH2D("display_empty_axis", "", 100, -100., 100., 100, -100., 100.);
    empty->SetTitle(Form("spill %d frame %d; x [mm]; y [mm]",
                         reader.spill_id(), reader.frame_index()));
    empty->SetStats(0);
    empty->Draw("AXIS");
    gPad->Update();
    return;
  }

  auto axis = new TH2D("display_axis", "", 100, -100., 100., 100, -100., 100.);
  axis->SetTitle(Form("spill %d frame %d Cherenkov hits; x [mm]; y [mm]",
                      reader.spill_id(), reader.frame_index()));
  axis->SetStats(0);
  axis->SetMinimum(color_min);
  axis->SetMaximum(color_max);
  axis->Draw("COLZ");
  axis->Reset();
  axis->GetZaxis()->SetTitle(use_reference ?
                             "time - reference [clock]" :
                             "time - min(time) [clock]");

  for (const auto &hit : drawable) {
    auto dt = use_reference ? hit.time - reference_time : hit.time - vmin;
    auto box = new TBox(hit.x - 0.5 * pixel_size, hit.y - 0.5 * pixel_size,
                        hit.x + 0.5 * pixel_size, hit.y + 0.5 * pixel_size);
    box->SetFillColor(color_index(dt, color_min, color_max));
    box->SetLineColor(kBlack);
    box->SetLineWidth(1);
    box->Draw("same");
  }

  auto label = new TText(-98., 94.,
                         use_reference ?
                         Form("hits=%zu  dt=[%.3f, %.3f]",
                              drawable.size(), vmin, vmax) :
                         Form("hits=%zu  time range=%.3f",
                              drawable.size(), vmax - vmin));
  label->SetTextSize(0.03);
  label->Draw("same");

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

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 900, 900);

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

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 900, 900);

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
                       double pixel_size)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 900, 900);

  while (reader.next_spill()) {
    while (reader.next_frame()) {
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

  auto canvas = new TCanvas("cDisplay", "Cherenkov frame display", 900, 900);

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

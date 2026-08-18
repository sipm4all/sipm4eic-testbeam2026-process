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

namespace {

struct map_cell_t {
  int device = 0;
  int fifo = 0;
  int column = 0;
  int pixel = 0;
  double x = 0.;
  double y = 0.;
  long long count = 0;
};

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

map_cell_t *
find_cell(std::vector<map_cell_t> &cells, const hit_t &hit)
{
  for (auto &cell : cells) {
    if (cell.device == hit.device &&
        cell.fifo == hit.fifo &&
        cell.column == hit.column &&
        cell.pixel == hit.pixel)
      return &cell;
  }

  map_cell_t cell;
  cell.device = hit.device;
  cell.fifo = hit.fifo;
  cell.column = hit.column;
  cell.pixel = hit.pixel;
  cell.x = hit.x;
  cell.y = hit.y;
  cells.push_back(cell);
  return &cells.back();
}

} // namespace

void
draw_map(const char *filename = "triggered.root",
         const char *outfilename = "map.root",
         double pixel_size = 3.2)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  std::vector<map_cell_t> cells;
  long long nframes = 0;
  long long nhits = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;
      for (const auto &hit : reader.cherenkov_hits()) {
        if (!finite_xy(hit))
          continue;
        auto cell = find_cell(cells, hit);
        ++cell->count;
        ++nhits;
      }
    }
  }

  auto fout = TFile::Open(outfilename, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  auto canvas = new TCanvas("cMap", "integrated Cherenkov hit map", 900, 900);
  canvas->SetRightMargin(0.16);
  gStyle->SetPalette(kViridis);

  if (cells.empty()) {
    auto empty = new TH2D("map_empty_axis", "", 100, -100., 100., 100, -100., 100.);
    empty->SetTitle("Cherenkov occupancy; x [mm]; y [mm]");
    empty->SetStats(0);
    empty->Draw("AXIS");
    empty->Write("hMapAxis");
    canvas->Write();
    fout->Close();
    std::cerr << "WARNING: no drawable hits found" << std::endl;
    return;
  }

  long long max_count = 0;

  for (const auto &cell : cells) {
    max_count = std::max(max_count, cell.count);
  }

  auto axis = new TH2D("hMapAxis", "", 100, -100., 100., 100, -100., 100.);
  axis->SetTitle("Cherenkov normalized occupancy; x [mm]; y [mm]");
  axis->SetStats(0);
  axis->SetMinimum(0.);
  axis->SetMaximum(1.);
  axis->Draw("COLZ");
  axis->Reset();
  axis->GetZaxis()->SetTitle("entries / max(entries)");

  for (const auto &cell : cells) {
    double occupancy = max_count > 0 ? double(cell.count) / double(max_count) : 0.;
    auto box = new TBox(cell.x - 0.5 * pixel_size, cell.y - 0.5 * pixel_size,
                        cell.x + 0.5 * pixel_size, cell.y + 0.5 * pixel_size);
    box->SetFillColor(color_index(occupancy, 0., 1.));
    box->SetLineColor(kBlack);
    box->SetLineWidth(1);
    box->Draw("same");
  }

  auto label = new TText(-98., 94.,
                         Form("frames=%lld  hits=%lld  channels=%zu",
                              nframes, nhits, cells.size()));
  label->SetTextSize(0.03);
  label->Draw("same");

  axis->Write();
  canvas->Write();
  fout->Close();

  std::cout << "frames processed: " << nframes << std::endl;
  std::cout << "drawable hits:     " << nhits << std::endl;
  std::cout << "channels:          " << cells.size() << std::endl;
  std::cout << "output:            " << outfilename << std::endl;
}

#include "../lib/frame_selection.h"

#include <TFile.h>
#include <TH2D.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace {

constexpr int n_theta_bins = 100;
constexpr int n_time_bins = 64;
constexpr double theta_min = 0.;
constexpr double theta_max = 0.1;
constexpr double time_min = -32.;
constexpr double time_max = 32.;
constexpr double resolution_theta = 0.002;
constexpr double resolution_time = 1.;
constexpr double gaussian_cut = 4.;

bool
passes_selections(const trigger_reader_t &reader,
                  const std::vector<selection_ptr_t> &selections)
{
  for (const auto &selection : selections) {
    if (!selection) {
      std::cerr << "ERROR: null selection" << std::endl;
      return false;
    }
    if (!selection->is_selected(reader))
      return false;
  }
  return true;
}

bool
valid_hit(const hit_t &hit)
{
  if (!std::isfinite(hit.x) || !std::isfinite(hit.y) ||
      !std::isfinite(hit.theta) || !std::isfinite(hit.time))
    return false;
  return hit.theta >= theta_min && hit.theta <= theta_max &&
         hit.time >= time_min && hit.time <= time_max;
}

} // namespace

void
irt_analysis(const char *filename = "ring.root",
                    const std::vector<selection_ptr_t> &selections = {},
                    const std::string &outfilename = "ring_theta_analysis.root",
                    const char *ring_name = "ring",
                    long long max_events = -1)
{
  trigger_reader_t reader;
  if (!reader.open(filename, ring_name ? ring_name : "ring"))
    return;
  if (!reader.has_rings()) {
    std::cerr << "ERROR: input must contain the ring tree" << std::endl;
    return;
  }

  std::unique_ptr<TFile> output(TFile::Open(outfilename.c_str(), "RECREATE"));
  if (!output || output->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  auto histogram = new TH2D(
      "hAverageCherenkovThetaVsNHits",
      "average IRT Cherenkov angle for one ring;Hough-selected hits;average #theta [rad]",
      128, 0., 128., 512, 0., 0.1);

  long long frames = 0;
  long long one_ring_frames = 0;
  long long filled = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      if (max_events >= 0 && frames >= max_events)
        break;
      ++frames;
      if (!passes_selections(reader, selections) || reader.rings().size() != 1)
        continue;

      std::vector<const hit_t *> hits;
      for (const auto &hit : reader.cherenkov_hits())
        if (valid_hit(hit))
          hits.push_back(&hit);

      double best_score = -1.;
      double best_theta = 0.;
      double best_time = 0.;
      for (int itheta = 0; itheta < n_theta_bins; ++itheta) {
        const double theta = theta_min + (itheta + 0.5) *
                             (theta_max - theta_min) / n_theta_bins;
        for (int itime = 0; itime < n_time_bins; ++itime) {
          const double time = time_min + (itime + 0.5) *
                              (time_max - time_min) / n_time_bins;
          double score = 0.;
          for (const auto *hit : hits) {
            const double theta_pull = (hit->theta - theta) / resolution_theta;
            const double time_pull = (hit->time - time) / resolution_time;
            score += std::exp(-0.5 * (theta_pull * theta_pull +
                                      time_pull * time_pull));
          }
          if (score > best_score) {
            best_score = score;
            best_theta = theta;
            best_time = time;
          }
        }
      }

      double sum = 0.;
      int nhits = 0;
      for (const auto *hit : hits) {
        if (std::abs(hit->theta - best_theta) > gaussian_cut * resolution_theta ||
            std::abs(hit->time - best_time) > gaussian_cut * resolution_time)
          continue;
        sum += hit->theta;
        ++nhits;
      }
      ++one_ring_frames;
      if (nhits == 0)
        continue;
      histogram->Fill(nhits, sum / nhits);
      ++filled;
    }
  }

  output->cd();
  histogram->Write();
  output->Close();

  std::cout << "frames processed:        " << frames << std::endl;
  std::cout << "one-ring frames:         " << one_ring_frames << std::endl;
  std::cout << "histogram fills:         " << filled << std::endl;
  std::cout << "output:                  " << outfilename << std::endl;
}

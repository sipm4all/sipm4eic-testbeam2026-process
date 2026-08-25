#include "../lib/frame_selection.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

constexpr double ring_analysis_time_to_ns = 3.125;
// Keep these association definitions aligned with display.C.
constexpr double ring_analysis_match_tolerance = 5.;
constexpr double ring_analysis_match_time_window = 2.;

namespace {

bool
passes_selections(const trigger_reader_t &reader,
                  const std::vector<selection_ptr_t> &selections)
{
  for (const auto &selection : selections) {
    if (!selection) {
      std::cerr << "ERROR: null ring-analysis selection" << std::endl;
      return false;
    }
    if (!selection->is_selected(reader))
      return false;
  }
  return true;
}

}

void
ring_analysis(const char *filename = "ring.root",
              reference_ptr_t reference = nullptr,
              const std::vector<selection_ptr_t> &selections = {},
              const std::string &outfilename = "ring_analysis.root")
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;
  if (!reader.has_rings()) {
    std::cerr << "ERROR: input file does not contain the ring tree" << std::endl;
    return;
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename
              << std::endl;
    return;
  }

  auto hDeltaT = new TH1D(
      "hDeltaTRing", "hit time relative to ring time;#Deltat [ns];hits",
      1024, -32., 32.);
  auto hDistance = new TH1D(
      "hRingDistance", "signed radial residual from ring;residual [mm];hits",
      1024, -100., 100.);
  auto hDeltaTVsDistance = new TH2D(
      "hDeltaTRingVsDistance",
      "hit time relative to ring time vs signed radial residual;residual [mm];#Deltat [ns]",
      1024, -100., 100., 1024, -32., 32.);
  auto hDeltaTMatched = new TH1D(
      "hDeltaTRingMatched",
      "matched hit time relative to ring time;#Deltat [ns];hits",
      256, -8., 8.);
  auto hDistanceMatched = new TH1D(
      "hRingDistanceMatched", "matched signed radial residual;residual [mm];hits",
      512, -8., 8.);
  auto hDeltaTUnmatched = new TH1D(
      "hDeltaTRingUnmatched",
      "unmatched hit time relative to ring time;#Deltat [ns];hits",
      1024, -32., 32.);
  auto hDistanceUnmatched = new TH1D(
      "hRingDistanceUnmatched",
      "unmatched signed radial residual;residual [mm];hits",
      1024, -100., 100.);
  auto hDeltaTVsDistanceUnmatched = new TH2D(
      "hDeltaTRingVsDistanceUnmatched",
      "unmatched hit time vs signed radial residual;residual [mm];#Deltat [ns]",
      1024, -100., 100., 1024, -32., 32.);

  long long nframes = 0;
  long long nselected_frames = 0;
  long long nrings = 0;
  long long nhits = 0;
  long long nmatched_hits = 0;
  long long nunmatched_hits = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;
      if (!passes_selections(reader, selections))
        continue;

      if (reference && !reference->process(reader))
        continue;

      ++nselected_frames;
      for (const auto &ring : reader.rings()) {
        if (!std::isfinite(ring.time) || ring.radius <= 0.)
          continue;
        ++nrings;

        for (const auto &hit : reader.cherenkov_hits()) {
          if (!std::isfinite(hit.time) || !std::isfinite(hit.x) ||
              !std::isfinite(hit.y))
            continue;

          const double delta_t_native = hit.time - ring.time;
          const double delta_t_ns = delta_t_native * ring_analysis_time_to_ns;
          const double signed_distance =
              frame_selection_detail::signed_ring_residual(ring, hit);
          if (!std::isfinite(signed_distance))
            continue;
          const double distance = std::abs(signed_distance);

          hDeltaT->Fill(delta_t_ns);
          hDistance->Fill(signed_distance);
          hDeltaTVsDistance->Fill(signed_distance, delta_t_ns);
          ++nhits;

          const bool matched =
              distance <= ring_analysis_match_tolerance &&
              std::abs(delta_t_native) <= ring_analysis_match_time_window;
          if (matched) {
            hDeltaTMatched->Fill(delta_t_ns);
            hDistanceMatched->Fill(signed_distance);
            ++nmatched_hits;
          } else {
            hDeltaTUnmatched->Fill(delta_t_ns);
            hDistanceUnmatched->Fill(signed_distance);
            hDeltaTVsDistanceUnmatched->Fill(signed_distance, delta_t_ns);
            ++nunmatched_hits;
          }
        }
      }
    }
  }

  fout->cd();
  hDeltaT->Write();
  hDistance->Write();
  hDeltaTVsDistance->Write();
  hDeltaTMatched->Write();
  hDistanceMatched->Write();
  hDeltaTUnmatched->Write();
  hDistanceUnmatched->Write();
  hDeltaTVsDistanceUnmatched->Write();
  fout->Close();

  std::cout << "frames processed:        " << nframes << std::endl;
  std::cout << "frames selected:         " << nselected_frames << std::endl;
  std::cout << "rings analyzed:          " << nrings << std::endl;
  std::cout << "hit/ring pairs:          " << nhits << std::endl;
  std::cout << "matched hit/ring pairs:  " << nmatched_hits << std::endl;
  std::cout << "unmatched hit/ring pairs:" << nunmatched_hits << std::endl;
  std::cout << "output:                  " << outfilename << std::endl;
}

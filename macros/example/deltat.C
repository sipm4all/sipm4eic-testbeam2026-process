#include "../lib/frame_selection.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

constexpr int deltat_nbins = 2048;
constexpr double deltat_min = -32.;
constexpr double deltat_max = 32.;
constexpr int deltat_nbins_narrow = 128;
constexpr double deltat_narrow_min = -2.;
constexpr double deltat_narrow_max = 2.;
constexpr int spill_nbins = 100;
constexpr double spill_min = 0.;
constexpr double spill_max = 100.;

struct delta_histograms_t {
  TH2D *delta = nullptr;
  TH2D *spill = nullptr;
  TH2D *tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  TH1D *ring_radius = nullptr;
  TH1D *ring_inliers = nullptr;
  TH1D *ring_residual = nullptr;
  TH2D *ring_center = nullptr;
  std::map<std::pair<int, int>, TH2D *> spill_by_fifo;

  delta_histograms_t(bool with_spill_tdc,
                     const ring_target_t *ring_target = nullptr)
  {
    constexpr int nchannels = 9 * 256;
    delta = new TH2D("hDeltaT", "", nchannels, 0., nchannels,
                     deltat_nbins, deltat_min, deltat_max);
    if (with_spill_tdc) {
      spill = new TH2D("hDeltaT_spill", "", spill_nbins, spill_min, spill_max,
                       deltat_nbins_narrow, deltat_narrow_min, deltat_narrow_max);
      for (int itdc = 0; itdc < 4; ++itdc)
        tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "",
                             256, 0., 256., deltat_nbins_narrow,
                             deltat_narrow_min, deltat_narrow_max);
    }
    if (ring_target) {
      ring_radius = new TH1D("hRingRadius", "stored ring radius;radius [mm];frames",
                             200, ring_target->min_radius, ring_target->max_radius);
      ring_inliers = new TH1D("hRingInliers", "stored ring inliers;inlier hits;frames",
                              256, -0.5, 255.5);
      const double residual_max = std::max(20., 4. * ring_target->tolerance);
      ring_residual = new TH1D("hRingResidual",
                               Form("stored-ring radial residual (selection <= %.3g mm);residual [mm];hits",
                                    ring_target->tolerance),
                               200, 0., residual_max);
      ring_center = new TH2D("hRingCenter",
                             "stored ring centre;x centre [mm];y centre [mm]",
                             200, -200., 200., 200, -200., 200.);
    }
  }

  void fill(const hit_t &target, double reference_time,
            int spill_id, int &nfills)
  {
    const double delta_t = target.time - reference_time;
    delta->Fill(frame_selection_detail::channel_number(target), delta_t);
    if (spill)
      spill->Fill(spill_id, delta_t);
    if (target.tdc >= 0 && target.tdc < 4 && tdc[target.tdc])
      tdc[target.tdc]->Fill(target.fine, delta_t);
    ++nfills;
  }

  void fill_all_fifos(const trigger_reader_t &reader,
                      double reference_time,
                      int spill_id,
                      int &nfills)
  {
    auto fill_hit = [&](const hit_t &hit) {
      if (hit.type != 1)
        return;

      const auto key = std::make_pair(hit.device, hit.fifo);
      auto found = spill_by_fifo.find(key);
      if (found == spill_by_fifo.end()) {
        auto name = Form("hDeltaT_spill_device%d_fifo%d", hit.device, hit.fifo);
        auto title = Form("spill delta-t;spill;delta-t [native units] (%d,%d)",
                          hit.device, hit.fifo);
        auto *hist = new TH2D(name, title, spill_nbins, spill_min, spill_max,
                              deltat_nbins_narrow, deltat_narrow_min,
                              deltat_narrow_max);
        found = spill_by_fifo.emplace(key, hist).first;
      }
      found->second->Fill(spill_id, hit.time - reference_time);
      ++nfills;
    };

    for (const auto &hit : reader.timing_hits())
      fill_hit(hit);
    for (const auto &hit : reader.cherenkov_hits())
      fill_hit(hit);
  }

  void write(double normalization, double fifo_normalization)
  {
    auto write_scaled = [normalization](TH1 *hist) {
      if (!hist)
        return;
      hist->Sumw2();
      if (normalization > 0.)
        hist->Scale(1. / normalization);
      hist->Write();
    };
    write_scaled(delta);
    write_scaled(spill);
    for (auto hist : tdc)
      write_scaled(hist);
    write_scaled(ring_radius);
    write_scaled(ring_inliers);
    write_scaled(ring_residual);
    write_scaled(ring_center);
    for (const auto &entry : spill_by_fifo)
      write_scaled(entry.second);
  }
};

std::vector<const hit_t *>
all_hits(const trigger_reader_t &reader)
{
  std::vector<const hit_t *> result;
  for (const auto &hit : reader.trigger_hits()) result.push_back(&hit);
  for (const auto &hit : reader.timing_hits()) result.push_back(&hit);
  for (const auto &hit : reader.cherenkov_hits()) result.push_back(&hit);
  return result;
}

bool
process_targets(const trigger_reader_t &reader,
                const std::vector<target_ptr_t> &targets)
{
  if (targets.empty()) {
    std::cerr << " --- no target selectors were supplied" << std::endl;
    return false;
  }

  for (const auto &target : targets) {
    if (!target) {
      std::cerr << " --- null target selector" << std::endl;
      return false;
    }
    if (!target->process(reader))
      return false;
  }
  return true;
}

bool
deltat_passes_selections(const trigger_reader_t &reader,
                         const std::vector<selection_ptr_t> &selections)
{
  for (const auto &selection : selections) {
    if (!selection) {
      std::cerr << " --- null event selection" << std::endl;
      return false;
    }
    if (!selection->is_selected(reader))
      return false;
  }
  return true;
}

bool
matches_all_targets(const hit_t &hit,
                    const std::vector<target_ptr_t> &targets)
{
  for (const auto &target : targets) {
    if (!target->matches(hit))
      return false;
  }
  return true;
}

const ring_target_t *
find_ring_target(const std::vector<target_ptr_t> &targets)
{
  for (const auto &target : targets) {
    if (const auto ring_target = dynamic_cast<const ring_target_t *>(target.get()))
      return ring_target;
  }
  return nullptr;
}

void
fill_ring_diagnostics(delta_histograms_t &histograms,
                      const ring_target_t *ring_target,
                      const trigger_reader_t &reader)
{
  if (!ring_target)
    return;

  for (const auto &ring : ring_target->accepted_rings) {
    histograms.ring_radius->Fill(ring.radius);
    histograms.ring_inliers->Fill(ring.ninliers);
    histograms.ring_center->Fill(ring.x0, ring.y0);
  }

  for (const auto &hit : reader.cherenkov_hits()) {
    double residual = std::numeric_limits<double>::infinity();
    for (const auto &ring : ring_target->accepted_rings)
      residual = std::min(residual,
                          frame_selection_detail::ring_residual(ring, hit));
    if (std::isfinite(residual))
      histograms.ring_residual->Fill(residual);
  }
}

void
deltat(const std::string &filename,
       const std::vector<target_ptr_t> &targets,
       reference_ptr_t reference,
       const std::vector<selection_ptr_t> &selections = {},
       const std::string &outfilename = "deltat.root")
{
  if (!reference) {
    std::cerr << " --- no reference selector was supplied" << std::endl;
    return;
  }

  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  const ring_target_t *ring_target = find_ring_target(targets);
  delta_histograms_t histograms(true, ring_target);
  int nframes_used = 0;
  int nframes_with_reference = 0;
  int nfills = 0;
  int nfifo_fills = 0;
  int nring_frames = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      if (!deltat_passes_selections(reader, selections))
        continue;
      if (!process_targets(reader, targets))
        continue;
      if (!reference->process(reader))
        continue;

      if (ring_target) {
        ++nring_frames;
        fill_ring_diagnostics(histograms, ring_target, reader);
      }

      const double reference_time = reference->reference_time();
      if (!std::isfinite(reference_time))
        continue;
      ++nframes_with_reference;
      histograms.fill_all_fifos(reader, reference_time,
                                reader.spill_id(), nfifo_fills);

      bool frame_has_target = false;
      for (const auto *hit : all_hits(reader)) {
        if (!matches_all_targets(*hit, targets))
          continue;
        frame_has_target = true;
        histograms.fill(*hit, reference_time, reader.spill_id(), nfills);
      }
      if (frame_has_target)
        ++nframes_used;
    }
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    return;
  }
  if (nframes_used == 0)
    std::cerr << " --- no frames with matching targets and reference found"
              << std::endl;
  std::cout << " --- frames with target/reference hits: " << nframes_used
            << std::endl;
  std::cout << " --- frames with valid reference: " << nframes_with_reference
            << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;
  std::cout << " --- all-FIFO histogram fills: " << nfifo_fills << std::endl;
  std::cout << " --- all-FIFO histograms: " << histograms.spill_by_fifo.size()
            << std::endl;
  if (ring_target)
    std::cout << " --- frames with accepted stored rings: " << nring_frames
              << std::endl;
  histograms.write(nframes_used, nframes_with_reference);
  fout->Close();
}

void
deltat(const std::string &filename,
       target_ptr_t target,
       reference_ptr_t reference,
       const std::vector<selection_ptr_t> &selections = {},
       const std::string &outfilename = "deltat.root")
{
  deltat(filename, std::vector<target_ptr_t>{target}, reference,
         selections, outfilename);
}

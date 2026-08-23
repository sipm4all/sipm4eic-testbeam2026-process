#include "../lib/trigger_reader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

constexpr int deltat_nbins = 2048;
constexpr double deltat_min = -32.;
constexpr double deltat_max = 32.;
constexpr int spill_nbins = 100;
constexpr double spill_min = 0.;
constexpr double spill_max = 100.;

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

struct timing_reference_t {
  std::string branch;

  timing_reference_t(const std::string &branch_ = "T")
    : branch(branch_)
  {
  }
};

struct ring_selection_t {
  int min_inliers;
  int max_inliers;
  int iterations;
  double tolerance;
  double min_radius;
  double max_radius;

  ring_selection_t(int min_inliers_ = 8,
                   int max_inliers_ = std::numeric_limits<int>::max(),
                   int iterations_ = 256,
                   double tolerance_ = 3.5,
                   double min_radius_ = 1.,
                   double max_radius_ = 200.)
    : min_inliers(min_inliers_), max_inliers(max_inliers_), iterations(iterations_),
      tolerance(tolerance_), min_radius(min_radius_), max_radius(max_radius_)
  {
  }
};

bool
match_field(int value, int requested)
{
  return requested < 0 || value == requested;
}

int
channel_number(const hit_t &hit)
{
  // The old format represented trigger words with a synthetic fifo=32.
  // Trigger hits have no fifo in the new tree, but they are not normally used
  // as histogram targets; retain that convention for compatibility.
  const int fifo = hit.type == 9 ? 32 : hit.fifo;
  const int chip = fifo / 4;
  const int local = hit.pixel + 4 * hit.column + 32 * chip;
  return local + 256 * (hit.device - 192);
}

bool
matches(const hit_t &hit, field_selector_t selector)
{
  return match_field(hit.type, selector.type) &&
         match_field(hit.device, selector.device) &&
         match_field(hit.fifo, selector.fifo) &&
         match_field(hit.column, selector.column) &&
         match_field(hit.pixel, selector.pixel);
}

bool
matches(const hit_t &hit, channel_selector_t selector)
{
  if (selector.type == 9) {
    // A type-9 channel selector is {type, trigger-device}.
    return hit.type == 9 && match_field(hit.device, selector.channel);
  }
  return hit.type == selector.type &&
         match_field(channel_number(hit), selector.channel);
}

struct fit_ring_t {
  double x = 0.;
  double y = 0.;
  double radius = 0.;
  int inliers = 0;
  double residual_sum = std::numeric_limits<double>::infinity();
  bool valid = false;

  double residual(double px, double py) const
  {
    return std::abs(std::hypot(px - x, py - y) - radius);
  }

  bool contains(double px, double py, double tolerance) const
  {
    return valid && residual(px, py) <= tolerance;
  }
};

struct ring_point_t {
  double x;
  double y;
};

bool
circle_from_three(const ring_point_t &a,
                  const ring_point_t &b,
                  const ring_point_t &c,
                  fit_ring_t &ring)
{
  const double determinant = 2. * (a.x * (b.y - c.y) +
                                   b.x * (c.y - a.y) +
                                   c.x * (a.y - b.y));
  if (std::abs(determinant) < 1.e-9)
    return false;

  const double aa = a.x * a.x + a.y * a.y;
  const double bb = b.x * b.x + b.y * b.y;
  const double cc = c.x * c.x + c.y * c.y;
  ring.x = (aa * (b.y - c.y) + bb * (c.y - a.y) + cc * (a.y - b.y)) / determinant;
  ring.y = (aa * (c.x - b.x) + bb * (a.x - c.x) + cc * (b.x - a.x)) / determinant;
  ring.radius = std::hypot(a.x - ring.x, a.y - ring.y);
  return std::isfinite(ring.x) && std::isfinite(ring.y) && std::isfinite(ring.radius);
}

fit_ring_t
fit_ring(const std::vector<hit_t> &hits,
         int spill,
         int frame,
         const ring_selection_t &selection)
{
  fit_ring_t best;
  std::vector<ring_point_t> points;
  for (const auto &hit : hits) {
    if (std::isfinite(hit.x) && std::isfinite(hit.y))
      points.push_back({hit.x, hit.y});
  }
  if (points.size() < 3)
    return best;

  std::mt19937 generator(static_cast<unsigned>(0x9e3779b9u ^
                                                static_cast<unsigned>(spill * 65537 + frame)));
  std::uniform_int_distribution<int> pick(0, static_cast<int>(points.size() - 1));
  for (int iteration = 0; iteration < selection.iterations; ++iteration) {
    const int ia = pick(generator);
    const int ib = pick(generator);
    const int ic = pick(generator);
    if (ia == ib || ia == ic || ib == ic)
      continue;

    fit_ring_t candidate;
    if (!circle_from_three(points[ia], points[ib], points[ic], candidate))
      continue;
    if (candidate.radius < selection.min_radius ||
        candidate.radius > selection.max_radius)
      continue;

    int inliers = 0;
    double residual = 0.;
    for (const auto &point : points) {
      const double distance = candidate.residual(point.x, point.y);
      if (distance <= selection.tolerance) {
        ++inliers;
        residual += distance;
      }
    }
    if (inliers > best.inliers ||
        (inliers == best.inliers && residual < best.residual_sum)) {
      candidate.inliers = inliers;
      candidate.residual_sum = residual;
      best = candidate;
    }
  }

  best.valid = best.inliers >= selection.min_inliers &&
               best.inliers <= selection.max_inliers;
  return best;
}

struct delta_histograms_t {
  TH2D *delta = nullptr;
  TH2D *spill = nullptr;
  TH2D *tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  TH1D *ring_radius = nullptr;
  TH1D *ring_inliers = nullptr;
  TH1D *ring_residual = nullptr;
  TH2D *ring_center = nullptr;

  delta_histograms_t(bool with_spill_tdc,
                     const ring_selection_t *ring_selection = nullptr)
  {
    constexpr int nchannels = 9 * 256;
    delta = new TH2D("hDeltaT", "", nchannels, 0., nchannels,
                     deltat_nbins, deltat_min, deltat_max);
    if (with_spill_tdc) {
      spill = new TH2D("hDeltaT_spill", "", spill_nbins, spill_min, spill_max,
                       deltat_nbins, deltat_min, deltat_max);
      for (int itdc = 0; itdc < 4; ++itdc)
        tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "",
                             256, 0., 256., deltat_nbins,
                             deltat_min, deltat_max);
    }
    if (ring_selection) {
      ring_radius = new TH1D("hRingRadius", "RANSAC ring radius;radius [mm];frames",
                             200, ring_selection->min_radius, ring_selection->max_radius);
      ring_inliers = new TH1D("hRingInliers", "RANSAC ring inliers;inlier hits;frames",
                              256, -0.5, 255.5);
      const double residual_max = std::max(20., 4. * ring_selection->tolerance);
      ring_residual = new TH1D("hRingResidual",
                               Form("RANSAC radial residual (selection <= %.3g mm);residual [mm];hits",
                                    ring_selection->tolerance),
                               200, 0., residual_max);
      ring_center = new TH2D("hRingCenter",
                             "RANSAC ring centre;x centre [mm];y centre [mm]",
                             200, -200., 200., 200, -200., 200.);
    }
  }

  void fill(const hit_t &target, double reference_time,
            int spill_id, int &nfills)
  {
    const double delta_t = target.time - reference_time;
    delta->Fill(channel_number(target), delta_t);
    if (spill)
      spill->Fill(spill_id, delta_t);
    if (target.tdc >= 0 && target.tdc < 4 && tdc[target.tdc])
      tdc[target.tdc]->Fill(target.fine, delta_t);
    ++nfills;
  }

  void write(double normalization)
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

double
timing_reference_value(const trigger_reader_t &reader,
                       const timing_reference_t &reference)
{
  if (reference.branch == "T0") return reader.T0();
  if (reference.branch == "T1") return reader.T1();
  if (reference.branch == "T") return reader.T();
  return std::numeric_limits<double>::quiet_NaN();
}

void
deltat(const std::string filename,
       field_selector_t reference,
       const std::string outfilename = "deltat.root")
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  delta_histograms_t histograms(false);
  int ntriggers = 0;
  int nfills = 0;
  while (reader.next_spill()) {
    while (reader.next_frame()) {
      const auto hits = all_hits(reader);
      const hit_t *reference_hit = nullptr;
      for (const auto *hit : hits) {
        if (matches(*hit, reference)) {
          reference_hit = hit;
          break;
        }
      }
      if (!reference_hit)
        continue;
      ++ntriggers;
      for (const auto *hit : hits) {
        if (hit == reference_hit)
          continue;
        histograms.fill(*hit, reference_hit->time, reader.spill_id(), nfills);
      }
    }
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    return;
  }
  if (ntriggers == 0)
    std::cerr << " --- no matching trigger hit found inside the stored frames" << std::endl;
  std::cout << " --- triggers found: " << ntriggers << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;
  histograms.write(ntriggers);
  fout->Close();
}

void
deltat_timing_reference_impl(const std::string filename,
                             channel_selector_t target_selector,
                             timing_reference_t timing_reference,
                             const ring_selection_t *ring_selection,
                             const std::string outfilename)
{
  if (timing_reference.branch != "T" && timing_reference.branch != "T0" &&
      timing_reference.branch != "T1") {
    std::cerr << " --- unsupported timing reference branch: " << timing_reference.branch
              << " (valid choices are T, T0, T1)" << std::endl;
    return;
  }

  trigger_reader_t reader;
  if (!reader.open(filename))
    return;
  if (!reader.has_timing()) {
    std::cerr << " --- input file does not contain timing estimator branches" << std::endl;
    return;
  }

  delta_histograms_t histograms(true, ring_selection);
  int nframes_used = 0;
  int nfills = 0;
  int nring_frames = 0;
  long long nring_target_candidates = 0;
  long long nring_target_selected = 0;
  long long nring_target_rejected = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      if (!reader.timing_valid())
        continue;
      const double reference_time = timing_reference_value(reader, timing_reference);
      if (!std::isfinite(reference_time))
        continue;

      fit_ring_t ring;
      if (ring_selection) {
        ring = fit_ring(reader.cherenkov_hits(), reader.spill_id(),
                        reader.frame_index(), *ring_selection);
        if (!ring.valid)
          continue;
        ++nring_frames;
        histograms.ring_radius->Fill(ring.radius);
        histograms.ring_inliers->Fill(ring.inliers);
        histograms.ring_center->Fill(ring.x, ring.y);
        for (const auto &hit : reader.cherenkov_hits()) {
          if (std::isfinite(hit.x) && std::isfinite(hit.y))
            histograms.ring_residual->Fill(ring.residual(hit.x, hit.y));
        }
      }

      std::vector<const hit_t *> targets;
      for (const auto *hit : all_hits(reader)) {
        if (!matches(*hit, target_selector))
          continue;
        if (ring_selection) {
          ++nring_target_candidates;
          const bool selected = hit->type == 1 &&
                                std::isfinite(hit->x) && std::isfinite(hit->y) &&
                                ring.contains(hit->x, hit->y, ring_selection->tolerance);
          if (!selected) {
            ++nring_target_rejected;
            continue;
          }
          ++nring_target_selected;
        }
        targets.push_back(hit);
      }
      if (targets.empty())
        continue;

      ++nframes_used;
      for (const auto *target : targets)
        histograms.fill(*target, reference_time, reader.spill_id(), nfills);
    }
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    return;
  }
  if (nframes_used == 0)
    std::cerr << " --- no valid timing-reference frames with target hits found" << std::endl;
  std::cout << " --- timing reference: " << timing_reference.branch << std::endl;
  std::cout << " --- frames with target/timing-reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;
  if (ring_selection) {
    std::cout << " --- frames with accepted Cherenkov rings: " << nring_frames << std::endl;
    std::cout << " --- ring target candidates: " << nring_target_candidates << std::endl;
    std::cout << " --- ring target hits selected: " << nring_target_selected << std::endl;
    std::cout << " --- ring target hits rejected: " << nring_target_rejected << std::endl;
  }
  histograms.write(nframes_used);
  fout->Close();
}

void
deltat(const std::string filename,
       channel_selector_t target_selector,
       timing_reference_t timing_reference,
       const std::string outfilename = "deltat.root")
{
  deltat_timing_reference_impl(filename, target_selector, timing_reference,
                                nullptr, outfilename);
}

void
deltat(const std::string filename,
       channel_selector_t target_selector,
       timing_reference_t timing_reference,
       ring_selection_t ring_selection,
       const std::string outfilename = "deltat.root")
{
  deltat_timing_reference_impl(filename, target_selector, timing_reference,
                                &ring_selection, outfilename);
}

void
deltat(const std::string filename,
       field_selector_t target_selector,
       field_selector_t reference_selector,
       const std::string outfilename = "deltat.root")
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  delta_histograms_t histograms(true);
  int nframes_used = 0;
  int nfills = 0;
  while (reader.next_spill()) {
    while (reader.next_frame()) {
      std::vector<const hit_t *> targets;
      std::vector<const hit_t *> references;
      for (const auto *hit : all_hits(reader)) {
        if (matches(*hit, target_selector)) targets.push_back(hit);
        if (matches(*hit, reference_selector)) references.push_back(hit);
      }
      if (targets.empty() || references.empty())
        continue;
      ++nframes_used;
      for (const auto *target : targets)
        for (const auto *reference : references)
          if (target != reference)
            histograms.fill(*target, reference->time, reader.spill_id(), nfills);
    }
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    return;
  }
  if (nframes_used == 0)
    std::cerr << " --- no frames with both target and reference hits found" << std::endl;
  std::cout << " --- frames with target/reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;
  histograms.write(nframes_used);
  fout->Close();
}

void
deltat(const std::string filename,
       channel_selector_t target_selector,
       channel_selector_t reference_selector,
       const std::string outfilename = "deltat.root")
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  delta_histograms_t histograms(true);
  int nframes_used = 0;
  int nfills = 0;
  while (reader.next_spill()) {
    while (reader.next_frame()) {
      std::vector<const hit_t *> targets;
      std::vector<const hit_t *> references;
      for (const auto *hit : all_hits(reader)) {
        if (matches(*hit, target_selector)) targets.push_back(hit);
        if (matches(*hit, reference_selector)) references.push_back(hit);
      }
      if (targets.empty() || references.empty())
        continue;
      ++nframes_used;
      for (const auto *target : targets)
        for (const auto *reference : references)
          if (target != reference)
            histograms.fill(*target, reference->time, reader.spill_id(), nfills);
    }
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    return;
  }
  if (nframes_used == 0)
    std::cerr << " --- no frames with both target and reference hits found" << std::endl;
  std::cout << " --- frames with target/reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;
  histograms.write(nframes_used);
  fout->Close();
}

void
deltat(const std::string filename,
       int trigger_type,
       int trigger_device,
       int trigger_fifo,
       int trigger_column,
       int trigger_pixel,
       const std::string outfilename = "deltat.root")
{
  deltat(filename,
         field_selector_t{trigger_type, trigger_device, trigger_fifo,
                          trigger_column, trigger_pixel},
         outfilename);
}

void
deltat(const std::string filename,
       int target_type,
       int target_device,
       int target_fifo,
       int target_column,
       int target_pixel,
       int reference_type,
       int reference_device,
       int reference_fifo,
       int reference_column,
       int reference_pixel,
       const std::string outfilename = "deltat.root")
{
  deltat(filename,
         field_selector_t{target_type, target_device, target_fifo,
                          target_column, target_pixel},
         field_selector_t{reference_type, reference_device, reference_fifo,
                          reference_column, reference_pixel},
         outfilename);
}

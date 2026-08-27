#include <TFile.h>
#include <TKey.h>
#include <TH1D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <random>
#include <string>
#include <vector>

#include "common.h"

namespace {

constexpr int maxrings = 8;

struct options_t {
  std::string input;
  std::string output;
  std::string branch_name = "ring";
  int nrings = 1;
  int min_inliers = 8;
  float max_shared_fraction = 0.5f;
  float xmin = -100.f;
  float xmax = 100.f;
  float ymin = -100.f;
  float ymax = 100.f;
  float rmin = 1.f;
  float rmax = 200.f;
  float tmin = -32.f;
  float tmax = 32.f;
  float tsigma = 15.f;
  Long64_t max_events = -1;
  int ransac_iterations = 128;
  float ransac_tolerance = 5.f;
  float ransac_time_window = 5.f;
  float hough_center_window = 5.f;
  float hough_radius_window = 5.f;
  float hough_time_window = 2.f;
  float hough_x_step = 1.f;
  float hough_y_step = 1.f;
  float hough_radius_step = 1.f;
  float hough_time_step = 1.f;
};

struct point_t {
  float x;
  float y;
  float time;
};

struct circle_t {
  float x = 0.f;
  float y = 0.f;
  float radius = 0.f;
  float time = 0.f;
  int inliers = 0;
  float residual_sum = std::numeric_limits<float>::infinity();
  float score = 0.f;
  float seed_x = 0.f;
  float seed_y = 0.f;
  float seed_radius = 0.f;
  float seed_time = 0.f;
  bool has_seed = false;
  bool valid = false;
};

bool
circle_from_three(const point_t &a, const point_t &b, const point_t &c,
                  circle_t &circle)
{
  const float determinant = 2.f * (a.x * (b.y - c.y) +
                                   b.x * (c.y - a.y) +
                                   c.x * (a.y - b.y));
  if (std::fabs(determinant) < 1.e-6f)
    return false;

  const float aa = a.x * a.x + a.y * a.y;
  const float bb = b.x * b.x + b.y * b.y;
  const float cc = c.x * c.x + c.y * c.y;
  circle.x = (aa * (b.y - c.y) + bb * (c.y - a.y) +
              cc * (a.y - b.y)) / determinant;
  circle.y = (aa * (c.x - b.x) + bb * (a.x - c.x) +
              cc * (b.x - a.x)) / determinant;
  circle.radius = std::hypot(a.x - circle.x, a.y - circle.y);
  return std::isfinite(circle.x) && std::isfinite(circle.y) &&
         std::isfinite(circle.radius);
}

float
circle_residual(const circle_t &circle, const point_t &point)
{
  return std::fabs(std::hypot(point.x - circle.x, point.y - circle.y) -
                   circle.radius);
}

float
mean_time(const std::vector<point_t> &points,
          const std::vector<int> &indices)
{
  if (indices.empty())
    return std::numeric_limits<float>::quiet_NaN();
  float sum = 0.f;
  for (int index : indices)
    sum += points[index].time;
  return sum / indices.size();
}

bool
fit_circle(const std::vector<point_t> &points,
           const std::vector<int> &indices,
           circle_t &circle)
{
  if (indices.size() < 3)
    return false;

  double normal[3][3] = {};
  double rhs[3] = {};
  for (int index : indices) {
    const double x = points[index].x;
    const double y = points[index].y;
    const double row[3] = {x, y, 1.};
    const double value = x * x + y * y;
    for (int i = 0; i < 3; ++i) {
      rhs[i] += row[i] * value;
      for (int j = 0; j < 3; ++j)
        normal[i][j] += row[i] * row[j];
    }
  }

  for (int column = 0; column < 3; ++column) {
    int pivot = column;
    for (int row = column + 1; row < 3; ++row) {
      if (std::abs(normal[row][column]) >
          std::abs(normal[pivot][column]))
        pivot = row;
    }
    if (std::abs(normal[pivot][column]) < 1.e-12)
      return false;
    if (pivot != column) {
      for (int j = column; j < 3; ++j)
        std::swap(normal[column][j], normal[pivot][j]);
      std::swap(rhs[column], rhs[pivot]);
    }
    for (int row = column + 1; row < 3; ++row) {
      const double factor = normal[row][column] / normal[column][column];
      for (int j = column; j < 3; ++j)
        normal[row][j] -= factor * normal[column][j];
      rhs[row] -= factor * rhs[column];
    }
  }

  double solution[3] = {};
  for (int row = 2; row >= 0; --row) {
    double value = rhs[row];
    for (int column = row + 1; column < 3; ++column)
      value -= normal[row][column] * solution[column];
    solution[row] = value / normal[row][row];
  }

  circle.x = 0.5f * solution[0];
  circle.y = 0.5f * solution[1];
  const double radius_squared = solution[2] + circle.x * circle.x +
                                circle.y * circle.y;
  if (radius_squared <= 0.)
    return false;
  circle.radius = std::sqrt(radius_squared);
  return std::isfinite(circle.x) && std::isfinite(circle.y) &&
         std::isfinite(circle.radius);
}

std::vector<circle_t>
ransac_seeds(const std::vector<point_t> &points,
            const std::vector<int> &indices,
            const options_t &opt,
            std::mt19937 &generator)
{
  std::vector<circle_t> hypotheses;
  if (indices.size() < 3 || opt.ransac_iterations <= 0)
    return hypotheses;

  std::vector<int> temporal_indices;
  std::vector<int> accepted;
  for (int iteration = 0; iteration < opt.ransac_iterations; ++iteration) {
    std::uniform_int_distribution<std::size_t> pick(0, indices.size() - 1);
    const int ia = indices[pick(generator)];
    temporal_indices.clear();
    for (int index : indices) {
      if (std::fabs(points[index].time - points[ia].time) <=
          opt.ransac_time_window)
        temporal_indices.push_back(index);
    }
    if (temporal_indices.size() < 3)
      continue;

    std::uniform_int_distribution<std::size_t> temporal_pick(
        0, temporal_indices.size() - 1);
    const int ib = temporal_indices[temporal_pick(generator)];
    const int ic = temporal_indices[temporal_pick(generator)];
    if (ia == ib || ia == ic || ib == ic)
      continue;

    circle_t hypothesis;
    if (!circle_from_three(points[ia], points[ib], points[ic], hypothesis) ||
        hypothesis.x < opt.xmin || hypothesis.x > opt.xmax ||
        hypothesis.y < opt.ymin || hypothesis.y > opt.ymax ||
        hypothesis.radius < opt.rmin || hypothesis.radius > opt.rmax)
      continue;

    auto collect = [&](const circle_t &model, float ring_time,
                       std::vector<int> &out, float &residual_sum) {
      out.clear();
      residual_sum = 0.f;
      for (int index : indices) {
        const float residual = circle_residual(model, points[index]);
        if (residual <= opt.ransac_tolerance &&
            std::fabs(points[index].time - ring_time) <=
                opt.ransac_time_window) {
          out.push_back(index);
          residual_sum += residual;
        }
      }
    };

    float residual_sum = 0.f;
    hypothesis.time = points[ia].time;
    collect(hypothesis, hypothesis.time, accepted, residual_sum);
    hypothesis.time = mean_time(points, accepted);
    collect(hypothesis, hypothesis.time, accepted, residual_sum);
    if (!fit_circle(points, accepted, hypothesis))
      continue;
    for (int refit = 0; refit < 2; ++refit) {
      hypothesis.time = mean_time(points, accepted);
      collect(hypothesis, hypothesis.time, accepted, residual_sum);
      if (!fit_circle(points, accepted, hypothesis)) {
        accepted.clear();
        break;
      }
    }
    if (accepted.size() < 3 || !std::isfinite(hypothesis.time) ||
        hypothesis.x < opt.xmin || hypothesis.x > opt.xmax ||
        hypothesis.y < opt.ymin || hypothesis.y > opt.ymax ||
        hypothesis.radius < opt.rmin || hypothesis.radius > opt.rmax)
      continue;

    hypothesis.time = mean_time(points, accepted);
    collect(hypothesis, hypothesis.time, accepted, residual_sum);
    hypothesis.inliers = static_cast<int>(accepted.size());
    hypothesis.residual_sum = residual_sum;
    hypotheses.push_back(hypothesis);
  }

  std::vector<circle_t> seeds;
  seeds.reserve(opt.nrings);
  std::vector<int> remaining = indices;
  for (int iseed = 0;
       iseed < opt.nrings &&
       static_cast<int>(remaining.size()) >= opt.min_inliers;
       ++iseed) {
    int best_hypothesis = -1;
    int best_inliers = 0;
    float best_residual_sum = std::numeric_limits<float>::infinity();
    std::vector<int> best_inlier_indices;
    for (int ihypothesis = 0;
         ihypothesis < static_cast<int>(hypotheses.size()); ++ihypothesis) {
      const circle_t &hypothesis = hypotheses[ihypothesis];
      std::vector<int> inlier_indices;
      float residual_sum = 0.f;
      for (int index : remaining) {
        const float residual = circle_residual(hypothesis, points[index]);
        if (residual <= opt.ransac_tolerance &&
            std::fabs(points[index].time - hypothesis.time) <=
                opt.ransac_time_window) {
          inlier_indices.push_back(index);
          residual_sum += residual;
        }
      }
      const int ninliers = static_cast<int>(inlier_indices.size());
      if (ninliers < opt.min_inliers || ninliers < best_inliers ||
          (ninliers == best_inliers && residual_sum >= best_residual_sum))
        continue;
      best_hypothesis = ihypothesis;
      best_inliers = ninliers;
      best_residual_sum = residual_sum;
      best_inlier_indices.swap(inlier_indices);
    }
    if (best_hypothesis < 0)
      break;

    circle_t seed = hypotheses[best_hypothesis];
    seed.inliers = best_inliers;
    seed.residual_sum = best_residual_sum;
    seeds.push_back(seed);

    std::vector<int> next_remaining;
    next_remaining.reserve(remaining.size() - best_inlier_indices.size());
    for (int index : remaining) {
      if (std::find(best_inlier_indices.begin(),
                    best_inlier_indices.end(), index) ==
          best_inlier_indices.end())
        next_remaining.push_back(index);
    }
    remaining.swap(next_remaining);
  }
  return seeds;
}

void
evaluate_candidate(const std::vector<point_t> &points,
                   const std::vector<int> &indices,
                   const options_t &opt,
                   circle_t &circle,
                   std::vector<int> &inlier_indices)
{
  constexpr float spatial_resolution = 2.1f;
  inlier_indices.clear();
  float residual_sum = 0.f;
  for (int index : indices) {
    const float spatial = circle_residual(circle, points[index]);
    const float temporal = std::fabs(points[index].time - circle.time);
    if (spatial <= 4.f * spatial_resolution &&
        temporal <= 4.f * opt.tsigma) {
      inlier_indices.push_back(index);
      residual_sum += spatial;
    }
  }
  circle.inliers = static_cast<int>(inlier_indices.size());
  circle.residual_sum = residual_sum;
  circle.valid = circle.inliers >= opt.min_inliers;
}

int
count_shared_hits(const std::vector<int> &first,
                  const std::vector<int> &second)
{
  int shared = 0;
  std::size_t i = 0;
  std::size_t j = 0;
  while (i < first.size() && j < second.size()) {
    if (first[i] == second[j]) {
      ++shared;
      ++i;
      ++j;
    } else if (first[i] < second[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return shared;
}

bool
is_duplicate_candidate(const circle_t &first,
                       const circle_t &second,
                       const options_t &opt)
{
  return std::fabs(first.x - second.x) <= opt.hough_x_step &&
         std::fabs(first.y - second.y) <= opt.hough_y_step &&
         std::fabs(first.radius - second.radius) <= opt.hough_radius_step &&
         std::fabs(first.time - second.time) <= opt.hough_time_step;
}

struct local_boundary_t {
  bool x = false;
  bool y = false;
  bool radius = false;
  bool time = false;

  bool any() const { return x || y || radius || time; }
};

local_boundary_t
local_scan_boundary(const circle_t &candidate, const circle_t &seed,
                    float x_window, float y_window,
                    float radius_window, float time_window,
                    float x_step, float y_step, float radius_step,
                    float time_step)
{
  constexpr float epsilon = 1.e-6f;
  return {
      std::fabs(candidate.x - seed.x) >=
          x_window - 0.5f * x_step - epsilon,
      std::fabs(candidate.y - seed.y) >=
          y_window - 0.5f * y_step - epsilon,
      std::fabs(candidate.radius - seed.radius) >=
          radius_window - 0.5f * radius_step - epsilon,
      std::fabs(candidate.time - seed.time) >=
          time_window - 0.5f * time_step - epsilon};
}

bool
can_expand_window(float seed, float candidate, float window, float global_min,
                  float global_max, float step)
{
  constexpr float epsilon = 1.e-6f;
  const bool lower_boundary = candidate <= seed - window +
                                             0.5f * step + epsilon;
  const bool upper_boundary = candidate >= seed + window -
                                             0.5f * step - epsilon;
  return (lower_boundary && seed - window > global_min + epsilon) ||
         (upper_boundary && seed + window < global_max - epsilon);
}

int
number_of_bins(float half_width, float step)
{
  return std::max(1, static_cast<int>(std::ceil(2.f * half_width / step)) + 1);
}

void
set_local_hough_grid(data_t &data, const circle_t &seed, const options_t &opt,
                     float x_window, float y_window, float radius_window,
                     float time_window)
{
  data.min.x = seed.x - x_window -
               0.5f * opt.hough_x_step;
  data.max.x = seed.x + x_window +
               0.5f * opt.hough_x_step;
  data.min.y = seed.y - y_window -
               0.5f * opt.hough_y_step;
  data.max.y = seed.y + y_window +
               0.5f * opt.hough_y_step;
  data.min.r = seed.radius - radius_window -
               0.5f * opt.hough_radius_step;
  data.max.r = seed.radius + radius_window +
               0.5f * opt.hough_radius_step;
  data.min.t = seed.time - time_window -
               0.5f * opt.hough_time_step;
  data.max.t = seed.time + time_window +
               0.5f * opt.hough_time_step;
  data.bins = {
    number_of_bins(x_window, opt.hough_x_step),
    number_of_bins(y_window, opt.hough_y_step),
    number_of_bins(radius_window, opt.hough_radius_step),
    number_of_bins(time_window, opt.hough_time_step)
  };
  data.sigma.t = opt.tsigma;
}

bool
resize_host_hough_data(data_t &data, int required_size, int &capacity)
{
  if (required_size <= capacity)
    return true;

  const int required_grid_size = (required_size + 255) / 256;
  try {
    auto map_x = std::make_unique<float[]>(required_size);
    auto map_y = std::make_unique<float[]>(required_size);
    auto map_r = std::make_unique<float[]>(required_size);
    auto map_t = std::make_unique<float[]>(required_size);
    auto hough_h = std::make_unique<float[]>(required_size);
    auto hough_nh = std::make_unique<int[]>(required_size);
    auto hough_rh = std::make_unique<float[]>(required_grid_size);
    auto hough_rhi = std::make_unique<int[]>(required_grid_size);

    delete[] data.map.x;
    delete[] data.map.y;
    delete[] data.map.r;
    delete[] data.map.t;
    delete[] data.hough.h;
    delete[] data.hough.nh;
    delete[] data.hough.rh;
    delete[] data.hough.rhi;
    data.map.x = map_x.release();
    data.map.y = map_y.release();
    data.map.r = map_r.release();
    data.map.t = map_t.release();
    data.hough.h = hough_h.release();
    data.hough.nh = hough_nh.release();
    data.hough.rh = hough_rh.release();
    data.hough.rhi = hough_rhi.release();
    capacity = required_size;
    return true;
  } catch (const std::bad_alloc &) {
    return false;
  }
}

bool
has_branch(TTree *tree, const char *name)
{
  if (tree->GetBranch(name))
    return true;
  std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
  return false;
}

bool
valid_bins(int bins)
{
  return bins > 0 && bins <= 100000;
}

bool
copy_auxiliary_trees(TFile *input, TFile *output,
                     const std::string &ring_name)
{
  output->cd();
  TIter keys(input->GetListOfKeys());
  while (auto key = static_cast<TKey *>(keys())) {
    const char *name = key->GetName();
    if (std::string(name) == "frames" || std::string(name) == "trigger" ||
        std::string(name) == "timing" || std::string(name) == "cherenkov" ||
        std::string(name) == "spill_participation" ||
        std::string(name) == ring_name)
      continue;

    auto object = std::unique_ptr<TObject>(key->ReadObj());
    auto tree = dynamic_cast<TTree *>(object.get());
    if (!tree)
      continue;
    if (!tree->CloneTree(-1, "fast")) {
      std::cerr << "ERROR: failed to copy auxiliary tree '" << name << "'"
                << std::endl;
      return false;
    }
  }
  return true;
}

bool
run(const options_t &opt)
{
  auto input = TFile::Open(opt.input.c_str(), "READ");
  if (!input || input->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << opt.input << std::endl;
    return false;
  }

  auto frames_in = dynamic_cast<TTree *>(input->Get("frames"));
  auto cherenkov_in = dynamic_cast<TTree *>(input->Get("cherenkov"));
  if (!frames_in || !cherenkov_in) {
    std::cerr << "ERROR: input must contain 'frames' and 'cherenkov' trees"
              << std::endl;
    input->Close();
    return false;
  }
  if (!has_branch(cherenkov_in, "nhits") ||
      !has_branch(cherenkov_in, "x") || !has_branch(cherenkov_in, "y") ||
      !has_branch(cherenkov_in, "time")) {
    input->Close();
    return false;
  }

  const Long64_t total_entries = frames_in->GetEntries();
  if (cherenkov_in->GetEntries() != total_entries) {
    std::cerr << "ERROR: frames/cherenkov entry-count mismatch: frames="
              << total_entries << " cherenkov=" << cherenkov_in->GetEntries()
              << std::endl;
    input->Close();
    return false;
  }
  const Long64_t entries = opt.max_events < 0
                               ? total_entries
                               : std::min(opt.max_events, total_entries);

  TTreeReader reader(cherenkov_in);
  TTreeReaderValue<UShort_t> nhits(reader, "nhits");
  TTreeReaderArray<Float_t> x(reader, "x");
  TTreeReaderArray<Float_t> y(reader, "y");
  TTreeReaderArray<Float_t> time(reader, "time");

  auto output = TFile::Open(opt.output.c_str(), "RECREATE");
  if (!output || output->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << opt.output
              << std::endl;
    input->Close();
    return false;
  }
  output->cd();

  auto frames_out = frames_in->CloneTree(0);
  auto trigger_in = dynamic_cast<TTree *>(input->Get("trigger"));
  auto timing_in = dynamic_cast<TTree *>(input->Get("timing"));
  auto trigger_out = trigger_in ? trigger_in->CloneTree(0) : nullptr;
  auto timing_out = timing_in ? timing_in->CloneTree(0) : nullptr;
  auto cherenkov_out = cherenkov_in->CloneTree(0);
  if (!frames_out || !cherenkov_out || (trigger_in && !trigger_out) ||
      (timing_in && !timing_out)) {
    std::cerr << "ERROR: failed to create output frame trees" << std::endl;
    output->Close();
    input->Close();
    return false;
  }

  UChar_t nring = 0;
  Float_t ring_x0[maxrings] = {};
  Float_t ring_y0[maxrings] = {};
  Float_t ring_r[maxrings] = {};
  Float_t ring_e[maxrings] = {};
  Float_t ring_phi[maxrings] = {};
  Float_t ring_time[maxrings] = {};
  UShort_t ring_ninliers[maxrings] = {};
  auto ring_out = new TTree(opt.branch_name.c_str(),
                            "circle candidates, one entry per frame");
  ring_out->Branch("nring", &nring, "nring/b");
  ring_out->Branch("ring_x0", ring_x0, "ring_x0[nring]/F");
  ring_out->Branch("ring_y0", ring_y0, "ring_y0[nring]/F");
  ring_out->Branch("ring_r", ring_r, "ring_r[nring]/F");
  ring_out->Branch("ring_e", ring_e, "ring_e[nring]/F");
  ring_out->Branch("ring_phi", ring_phi, "ring_phi[nring]/F");
  ring_out->Branch("ring_time", ring_time, "ring_time[nring]/F");
  ring_out->Branch("ring_ninliers", ring_ninliers, "ring_ninliers[nring]/s");

  // Diagnostic distributions use final minus RANSAC-seed parameters. They
  // are histograms only; the persistent ring-tree schema is unchanged.
  auto hRansacDeltaX0 = new TH1D(
      "hRansacDeltaX0", "final x0 - RANSAC seed x0;#Delta x0;entries",
      200, -10., 10.);
  auto hRansacDeltaY0 = new TH1D(
      "hRansacDeltaY0", "final y0 - RANSAC seed y0;#Delta y0;entries",
      200, -10., 10.);
  auto hRansacDeltaRadius = new TH1D(
      "hRansacDeltaRadius",
      "final radius - RANSAC seed radius;#Delta R;entries", 200, -10., 10.);
  auto hRansacDeltaTime = new TH1D(
      "hRansacDeltaTime",
      "final ring time - RANSAC seed time;#Delta t;entries", 200, -5., 5.);

  // Allocate host-side maps for the initial scan and grow them if successive
  // boundary retries produce a larger per-seed grid.
  const int max_xbins = number_of_bins(2.f * opt.hough_center_window,
                                      opt.hough_x_step);
  const int max_ybins = number_of_bins(2.f * opt.hough_center_window,
                                      opt.hough_y_step);
  const int max_rbins = number_of_bins(2.f * opt.hough_radius_window,
                                      opt.hough_radius_step);
  const int max_tbins = number_of_bins(2.f * opt.hough_time_window,
                                       opt.hough_time_step);
  const int max_size = max_xbins * max_ybins * max_rbins * max_tbins;
  const int max_grid_size = (max_size + 255) / 256;
  data_t data{};
  data.map.x = new float[max_size];
  data.map.y = new float[max_size];
  data.map.r = new float[max_size];
  data.map.t = new float[max_size];
  data.hough.h = new float[max_size];
  data.hough.nh = new int[max_size];
  data.hough.rh = new float[max_grid_size];
  data.hough.rhi = new int[max_grid_size];
  int data_capacity = max_size;

  Long64_t processed = 0;
  Long64_t rings_found = 0;
  Long64_t ransac_seed_count = 0;
  Long64_t frames_with_seed = 0;
  Long64_t boundary_retries = 0;
  std::vector<float> points_x;
  std::vector<float> points_y;
  std::vector<float> points_t;
  points_x.reserve(1 << 16);
  points_y.reserve(1 << 16);
  points_t.reserve(1 << 16);
  std::vector<point_t> points;
  std::mt19937 generator(0x9e3779b9u);

  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (frames_in->GetEntry(entry) <= 0 || !reader.Next()) {
      std::cerr << "ERROR: failed to read frame entry " << entry << std::endl;
      hough_free();
      delete[] data.map.x;
      delete[] data.map.y;
      delete[] data.map.r;
      delete[] data.map.t;
      delete[] data.hough.h;
      delete[] data.hough.nh;
      delete[] data.hough.rh;
      delete[] data.hough.rhi;
      output->Close();
      input->Close();
      return false;
    }

    const int n = static_cast<int>(*nhits);
    if (x.GetSize() < n || y.GetSize() < n || time.GetSize() < n ||
        n > (1 << 16)) {
      std::cerr << "ERROR: invalid Cherenkov hit array at frame " << entry
                << " nhits=" << n << std::endl;
      hough_free();
      delete[] data.map.x;
      delete[] data.map.y;
      delete[] data.map.r;
      delete[] data.map.t;
      delete[] data.hough.h;
      delete[] data.hough.nh;
      delete[] data.hough.rh;
      delete[] data.hough.rhi;
      output->Close();
      input->Close();
      return false;
    }

    points_x.assign(x.begin(), x.begin() + n);
    points_y.assign(y.begin(), y.begin() + n);
    points_t.assign(time.begin(), time.begin() + n);
    data.points = {n, points_x.data(), points_y.data(), points_t.data()};
    points.resize(n);
    for (int hit = 0; hit < n; ++hit)
      points[hit] = {points_x[hit], points_y[hit], points_t[hit]};

    nring = 0;
    std::vector<int> active;
    active.reserve(n);
    for (int hit = 0; hit < n; ++hit) {
      if (data.points.x[hit] != 0.f && data.points.y[hit] != 0.f)
        active.push_back(hit);
    }
    const std::vector<circle_t> seeds =
        ransac_seeds(points, active, opt, generator);
    ransac_seed_count += seeds.size();
    if (!seeds.empty())
      ++frames_with_seed;

    // Each seed gets its own local accumulator and contributes exactly one
    // maximum. Candidate selection follows ring-finder-hough: sort maxima,
    // remove duplicate models, validate inliers, and reject excessive
    // sharing. Hits remain available to overlapping candidates.
    std::vector<circle_t> candidates;
    auto scan_seed = [&](const circle_t &seed, float x_window,
                         float y_window, float radius_window, float time_window,
                         circle_t &candidate) {
      candidate.has_seed = false;
      set_local_hough_grid(data, seed, opt, x_window, y_window,
                           radius_window,
                           time_window);
      const long long local_size_long =
          static_cast<long long>(data.bins.x) * data.bins.y * data.bins.r *
          data.bins.t;
      if (local_size_long <= 0 ||
          local_size_long > std::numeric_limits<int>::max() ||
          !resize_host_hough_data(data, static_cast<int>(local_size_long),
                                  data_capacity)) {
        std::cerr << "ERROR: adaptive Hough grid is too large" << std::endl;
        return false;
      }
      hough_init(data);
      hough_transform(data);
      const int local_size = data.bins.x * data.bins.y * data.bins.r *
                             data.bins.t;
      const int local_grid_size = (local_size + 255) / 256;
      const int block = static_cast<int>(std::distance(
          data.hough.rh,
          std::max_element(data.hough.rh,
                           data.hough.rh + local_grid_size)));
      const int maximum = data.hough.rhi[block];
      if (maximum < 0)
        return false;

      candidate.x = data.map.x[maximum];
      candidate.y = data.map.y[maximum];
      candidate.radius = data.map.r[maximum];
      candidate.time = data.map.t[maximum];
      candidate.score = data.hough.h[maximum];
      candidate.seed_x = seed.x;
      candidate.seed_y = seed.y;
      candidate.seed_radius = seed.radius;
      candidate.seed_time = seed.time;
      candidate.has_seed = true;
      return true;
    };
    for (const circle_t &seed : seeds) {
      circle_t candidate;
      float x_window = opt.hough_center_window;
      float y_window = opt.hough_center_window;
      float radius_window = opt.hough_radius_window;
      float time_window = opt.hough_time_window;
      for (;;) {
        if (!scan_seed(seed, x_window, y_window, radius_window, time_window,
                       candidate))
          break;
        const local_boundary_t boundary = local_scan_boundary(
            candidate, seed, x_window, y_window, radius_window, time_window,
            opt.hough_x_step, opt.hough_y_step, opt.hough_radius_step,
            opt.hough_time_step);
        bool expanded = false;
        if (boundary.x && can_expand_window(
                              seed.x, candidate.x, x_window, opt.xmin,
                              opt.xmax, opt.hough_x_step)) {
          x_window *= 2.f;
          expanded = true;
        }
        if (boundary.y && can_expand_window(
                              seed.y, candidate.y, y_window, opt.ymin,
                              opt.ymax, opt.hough_y_step)) {
          y_window *= 2.f;
          expanded = true;
        }
        if (boundary.radius && can_expand_window(
                                  seed.radius, candidate.radius,
                                  radius_window, opt.rmin, opt.rmax,
                                  opt.hough_radius_step)) {
          radius_window *= 2.f;
          expanded = true;
        }
        if (boundary.time && can_expand_window(
                                 seed.time, candidate.time, time_window,
                                 opt.tmin, opt.tmax, opt.hough_time_step)) {
          time_window *= 2.f;
          expanded = true;
        }
        if (!expanded)
          break;
        ++boundary_retries;
      }
      if (!candidate.has_seed)
        continue;
      candidates.push_back(candidate);
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const circle_t &first, const circle_t &second) {
                return first.score > second.score;
              });
    std::vector<circle_t> unique_candidates;
    unique_candidates.reserve(candidates.size());
    for (const circle_t &candidate : candidates) {
      bool duplicate = false;
      for (const circle_t &kept : unique_candidates) {
        if (is_duplicate_candidate(candidate, kept, opt)) {
          duplicate = true;
          break;
        }
      }
      if (!duplicate)
        unique_candidates.push_back(candidate);
    }

    std::vector<std::vector<int>> accepted_inliers;
    for (circle_t candidate : unique_candidates) {
      std::vector<int> inliers;
      evaluate_candidate(points, active, opt, candidate, inliers);
      if (!candidate.valid)
        continue;

      bool too_many_shared_hits = false;
      for (const auto &previous_inliers : accepted_inliers) {
        const int shared_hits = count_shared_hits(inliers, previous_inliers);
        const std::size_t smaller_ring =
            std::min(inliers.size(), previous_inliers.size());
        const int max_shared = static_cast<int>(std::floor(
            opt.max_shared_fraction * static_cast<float>(smaller_ring)));
        if (shared_hits > max_shared) {
          too_many_shared_hits = true;
          break;
        }
      }
      if (too_many_shared_hits)
        continue;

      ring_x0[nring] = candidate.x;
      ring_y0[nring] = candidate.y;
      ring_r[nring] = candidate.radius;
      ring_e[nring] = 0.f;
      ring_phi[nring] = 0.f;
      ring_time[nring] = candidate.time;
      ring_ninliers[nring] = static_cast<UShort_t>(candidate.inliers);
      if (candidate.has_seed) {
        hRansacDeltaX0->Fill(candidate.x - candidate.seed_x);
        hRansacDeltaY0->Fill(candidate.y - candidate.seed_y);
        hRansacDeltaRadius->Fill(candidate.radius - candidate.seed_radius);
        hRansacDeltaTime->Fill(candidate.time - candidate.seed_time);
      }
      ++nring;
      ++rings_found;
      accepted_inliers.push_back(inliers);
      if (nring == maxrings)
        break;
    }

    if (trigger_in)
      trigger_in->GetEntry(entry);
    if (timing_in)
      timing_in->GetEntry(entry);
    cherenkov_in->GetEntry(entry);
    frames_out->Fill();
    if (trigger_out)
      trigger_out->Fill();
    if (timing_out)
      timing_out->Fill();
    cherenkov_out->Fill();
    ring_out->Fill();
    ++processed;
  }

  hough_free();
  delete[] data.map.x;
  delete[] data.map.y;
  delete[] data.map.r;
  delete[] data.map.t;
  delete[] data.hough.h;
  delete[] data.hough.nh;
  delete[] data.hough.rh;
  delete[] data.hough.rhi;

  for (auto tree : {frames_out, trigger_out, timing_out, cherenkov_out,
                    ring_out}) {
    if (tree && tree->GetEntries() != processed) {
      std::cerr << "ERROR: output tree entry-count mismatch in "
                << tree->GetName() << std::endl;
      output->Close();
      input->Close();
      return false;
    }
  }

  if (auto participation = dynamic_cast<TTree *>(
          input->Get("spill_participation"))) {
    output->cd();
    if (!participation->CloneTree(-1, "fast")) {
      std::cerr << "ERROR: failed to copy spill_participation tree" << std::endl;
      output->Close();
      input->Close();
      return false;
    }
  }
  if (!copy_auxiliary_trees(input, output, opt.branch_name)) {
    output->Close();
    input->Close();
    return false;
  }

  output->cd();
  output->Write();
  output->Close();
  input->Close();
  std::cout << "frames processed: " << processed << std::endl
            << "RANSAC seeds:     " << ransac_seed_count << " ("
            << frames_with_seed << " ring iterations)" << std::endl
            << "boundary retries: " << boundary_retries << std::endl
            << "rings found:      " << rings_found << std::endl
            << "output:            " << opt.output << std::endl;
  return true;
}

}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;
  options_t opt;
  po::options_description description("options");
  description.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&opt.input)->required(),
     "triggered input ROOT file")
    ("output,o", po::value<std::string>(&opt.output)->required(),
     "output ROOT file")
    ("branch-name", po::value<std::string>(&opt.branch_name)
                         ->default_value(opt.branch_name),
     "output ring tree name")
    ("nrings", po::value<int>(&opt.nrings)->default_value(opt.nrings),
     "maximum number of rings per frame")
    ("min-inliers", po::value<int>(&opt.min_inliers)
                              ->default_value(opt.min_inliers),
     "minimum RANSAC inliers per seed")
    ("max-shared-fraction", po::value<float>(&opt.max_shared_fraction)
                                  ->default_value(opt.max_shared_fraction),
     "maximum fraction of the smaller ring's inliers shared by two rings")
    ("xmin", po::value<float>(&opt.xmin)->default_value(opt.xmin),
     "Hough x minimum")
    ("xmax", po::value<float>(&opt.xmax)->default_value(opt.xmax),
     "Hough x maximum")
    ("ymin", po::value<float>(&opt.ymin)->default_value(opt.ymin),
     "Hough y minimum")
    ("ymax", po::value<float>(&opt.ymax)->default_value(opt.ymax),
     "Hough y maximum")
    ("rmin", po::value<float>(&opt.rmin)->default_value(opt.rmin),
     "Hough radius minimum")
    ("rmax", po::value<float>(&opt.rmax)->default_value(opt.rmax),
     "Hough radius maximum")
    ("tmin", po::value<float>(&opt.tmin)->default_value(opt.tmin),
     "Hough time minimum")
    ("tmax", po::value<float>(&opt.tmax)->default_value(opt.tmax),
     "Hough time maximum")
    ("tsigma", po::value<float>(&opt.tsigma)->default_value(opt.tsigma),
     "Hough time resolution")
    ("ransac-iterations", po::value<int>(&opt.ransac_iterations)
                              ->default_value(opt.ransac_iterations),
     "RANSAC iterations per ring")
    ("ransac-tolerance", po::value<float>(&opt.ransac_tolerance)
                              ->default_value(opt.ransac_tolerance),
     "RANSAC spatial inlier tolerance")
    ("ransac-time-window", po::value<float>(&opt.ransac_time_window)
                                ->default_value(opt.ransac_time_window),
     "RANSAC temporal inlier window")
    ("hough-center-window", po::value<float>(&opt.hough_center_window)
                                  ->default_value(opt.hough_center_window),
     "Hough x/y half-window around the RANSAC seed")
    ("hough-radius-window", po::value<float>(&opt.hough_radius_window)
                                  ->default_value(opt.hough_radius_window),
     "Hough radius half-window around the RANSAC seed")
    ("hough-time-window", po::value<float>(&opt.hough_time_window)
                                ->default_value(opt.hough_time_window),
     "Hough time half-window around the RANSAC seed")
    ("hough-x-step", po::value<float>(&opt.hough_x_step)
                            ->default_value(opt.hough_x_step),
     "Hough x step")
    ("hough-y-step", po::value<float>(&opt.hough_y_step)
                            ->default_value(opt.hough_y_step),
     "Hough y step")
    ("hough-radius-step", po::value<float>(&opt.hough_radius_step)
                                 ->default_value(opt.hough_radius_step),
     "Hough radius step")
    ("hough-time-step", po::value<float>(&opt.hough_time_step)
                               ->default_value(opt.hough_time_step),
     "Hough time step")
    ("max-events", po::value<Long64_t>(&opt.max_events)
                        ->default_value(opt.max_events),
     "maximum frames to process; negative means all");

  try {
    po::variables_map variables;
    po::store(po::parse_command_line(argc, argv, description), variables);
    if (variables.count("help")) {
      std::cout << description << std::endl;
      return 0;
    }
    po::notify(variables);
    if (opt.input.empty() || opt.output.empty() || opt.branch_name.empty() ||
        opt.nrings < 1 || opt.nrings > maxrings || opt.min_inliers < 3 ||
        opt.max_shared_fraction < 0.f || opt.max_shared_fraction > 1.f ||
        opt.max_events < -1 ||
        opt.xmin > opt.xmax || opt.ymin > opt.ymax || opt.rmin >= opt.rmax ||
        opt.tmin >= opt.tmax || opt.tsigma <= 0.f ||
        opt.ransac_iterations <= 0 || opt.ransac_tolerance <= 0.f ||
        opt.ransac_time_window <= 0.f || opt.hough_center_window <= 0.f ||
        opt.hough_radius_window <= 0.f || opt.hough_time_window <= 0.f ||
        opt.hough_x_step <= 0.f || opt.hough_y_step <= 0.f ||
        opt.hough_radius_step <= 0.f || opt.hough_time_step <= 0.f ||
        !valid_bins(number_of_bins(opt.hough_center_window, opt.hough_x_step)) ||
        !valid_bins(number_of_bins(opt.hough_center_window, opt.hough_y_step)) ||
        !valid_bins(number_of_bins(opt.hough_radius_window,
                                   opt.hough_radius_step)) ||
        !valid_bins(number_of_bins(opt.hough_time_window,
                                   opt.hough_time_step)))
      throw std::runtime_error("invalid ring-finder-roberto option");
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << std::endl
              << description << std::endl;
    return 1;
  }

  return run(opt) ? 0 : 1;
}

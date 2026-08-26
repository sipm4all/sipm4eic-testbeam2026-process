#include <TFile.h>
#include <TKey.h>
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
  float hough_center_window = 10.f;
  float hough_radius_window = 10.f;
  float hough_time_window = 5.f;
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

bool
collect_ransac_inliers(const std::vector<point_t> &points,
                       const std::vector<int> &indices,
                       const circle_t &circle,
                       float ring_time,
                       float spatial_tolerance,
                       float time_window,
                       std::vector<int> &inliers,
                       float &residual_sum)
{
  inliers.clear();
  residual_sum = 0.f;
  for (int index : indices) {
    const float residual = circle_residual(circle, points[index]);
    if (residual > spatial_tolerance ||
        std::fabs(points[index].time - ring_time) > time_window)
      continue;
    inliers.push_back(index);
    residual_sum += residual;
  }
  return !inliers.empty();
}

bool
find_ransac_seed(const std::vector<point_t> &points,
                 const std::vector<int> &indices,
                 const options_t &opt,
                 std::mt19937 &generator,
                 circle_t &best)
{
  if (indices.size() < 3 || opt.ransac_iterations <= 0)
    return false;

  std::vector<int> temporal_indices;
  std::vector<int> inliers;
  int best_inliers = 0;
  float best_residual_sum = std::numeric_limits<float>::infinity();
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

    circle_t candidate;
    if (!circle_from_three(points[ia], points[ib], points[ic], candidate) ||
        candidate.x < opt.xmin || candidate.x > opt.xmax ||
        candidate.y < opt.ymin || candidate.y > opt.ymax ||
        candidate.radius < opt.rmin || candidate.radius > opt.rmax)
      continue;

    float residual_sum = 0.f;
    collect_ransac_inliers(points, indices, candidate, points[ia].time,
                           opt.ransac_tolerance, opt.ransac_time_window,
                           inliers, residual_sum);
    if (static_cast<int>(inliers.size()) < best_inliers ||
        (static_cast<int>(inliers.size()) == best_inliers &&
         residual_sum >= best_residual_sum))
      continue;

    best = candidate;
    best.time = points[ia].time;
    best.inliers = static_cast<int>(inliers.size());
    best.residual_sum = residual_sum;
    best.valid = true;
    best_inliers = best.inliers;
    best_residual_sum = residual_sum;
  }
  return best.valid;
}

int
number_of_bins(float half_width, float step)
{
  return std::max(1, static_cast<int>(std::ceil(2.f * half_width / step)) + 1);
}

void
set_local_hough_grid(data_t &data, const circle_t &seed, const options_t &opt)
{
  data.min.x = seed.x - opt.hough_center_window -
               0.5f * opt.hough_x_step;
  data.max.x = seed.x + opt.hough_center_window +
               0.5f * opt.hough_x_step;
  data.min.y = seed.y - opt.hough_center_window -
               0.5f * opt.hough_y_step;
  data.max.y = seed.y + opt.hough_center_window +
               0.5f * opt.hough_y_step;
  data.min.r = seed.radius - opt.hough_radius_window -
               0.5f * opt.hough_radius_step;
  data.max.r = seed.radius + opt.hough_radius_window +
               0.5f * opt.hough_radius_step;
  data.min.t = seed.time - opt.hough_time_window -
               0.5f * opt.hough_time_step;
  data.max.t = seed.time + opt.hough_time_window +
               0.5f * opt.hough_time_step;
  data.bins = {
    number_of_bins(opt.hough_center_window, opt.hough_x_step),
    number_of_bins(opt.hough_center_window, opt.hough_y_step),
    number_of_bins(opt.hough_radius_window, opt.hough_radius_step),
    number_of_bins(opt.hough_time_window, opt.hough_time_step)
  };
  data.sigma.t = opt.tsigma;
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

  const int xbins = number_of_bins(opt.hough_center_window,
                                   opt.hough_x_step);
  const int ybins = number_of_bins(opt.hough_center_window,
                                   opt.hough_y_step);
  const int rbins = number_of_bins(opt.hough_radius_window,
                                   opt.hough_radius_step);
  const int tbins = number_of_bins(opt.hough_time_window,
                                   opt.hough_time_step);
  const int size = xbins * ybins * rbins * tbins;
  const int grid_size = (size + 255) / 256;
  data_t data{};
  data.map.x = new float[size];
  data.map.y = new float[size];
  data.map.r = new float[size];
  data.map.t = new float[size];
  data.hough.h = new float[size];
  data.hough.nh = new int[size];
  data.hough.rh = new float[grid_size];
  data.hough.rhi = new int[grid_size];

  Long64_t processed = 0;
  Long64_t rings_found = 0;
  Long64_t ransac_seeds = 0;
  Long64_t frames_with_seed = 0;
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
    int used = 0;
    while (nring < opt.nrings && used < n) {
      std::vector<int> active;
      active.reserve(n - used);
      for (int hit = 0; hit < n; ++hit) {
        if (data.points.x[hit] != 0.f && data.points.y[hit] != 0.f)
          active.push_back(hit);
      }
      circle_t seed;
      if (!find_ransac_seed(points, active, opt, generator, seed))
        break;
      ++ransac_seeds;
      ++frames_with_seed;

      set_local_hough_grid(data, seed, opt);
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
        break;

      ring_x0[nring] = data.map.x[maximum];
      ring_y0[nring] = data.map.y[maximum];
      ring_r[nring] = data.map.r[maximum];
      ring_e[nring] = 0.f;
      ring_phi[nring] = 0.f;
      ring_time[nring] = data.map.t[maximum];
      ring_ninliers[nring] = static_cast<UShort_t>(data.hough.nh[maximum]);

      for (int hit = 0; hit < n; ++hit) {
        if (data.points.x[hit] == 0.f || data.points.y[hit] == 0.f)
          continue;
        const float dx = data.points.x[hit] - ring_x0[nring];
        const float dy = data.points.y[hit] - ring_y0[nring];
        const float radius = std::hypot(dx, dy);
        if (std::fabs(radius - ring_r[nring]) > 3.f * 2.1f)
          continue;
        if (std::fabs(data.points.t[hit] - ring_time[nring]) >
            3.f * opt.tsigma)
          continue;
        data.points.x[hit] = 0.f;
        data.points.y[hit] = 0.f;
        points[hit].x = 0.f;
        points[hit].y = 0.f;
        ++used;
      }
      ++nring;
      ++rings_found;
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
            << "RANSAC seeds:     " << ransac_seeds << " ("
            << frames_with_seed << " ring iterations)" << std::endl
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
        opt.nrings < 1 || opt.nrings > maxrings || opt.max_events < -1 ||
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

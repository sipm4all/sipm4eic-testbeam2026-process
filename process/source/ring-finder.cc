#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr int maxrings = 8;

struct point_t {
  double x;
  double y;
  double time;
};

struct circle_t {
  double x = 0.;
  double y = 0.;
  double radius = 0.;
  double eccentricity = 0.;
  double phi = 0.;
  double ring_time = 0.;
  int inliers = 0;
  double residual_sum = std::numeric_limits<double>::infinity();
  bool valid = false;
};

double ellipse_residual(const circle_t &ellipse, const point_t &point)
{
  double dx = point.x - ellipse.x;
  double dy = point.y - ellipse.y;
  double xp = std::cos(ellipse.phi) * dx + std::sin(ellipse.phi) * dy;
  double yp = -std::sin(ellipse.phi) * dx + std::cos(ellipse.phi) * dy;
  double b = ellipse.radius * std::sqrt(std::max(1.e-12,
                                                   1. - ellipse.eccentricity * ellipse.eccentricity));
  return std::abs(std::hypot(xp / ellipse.radius, yp / b) - 1.) * ellipse.radius;
}

void refine_ellipse(circle_t &ellipse, const std::vector<point_t> &points,
                    const std::vector<int> &indices, bool force_circle)
{
  if (force_circle || indices.size() < 5) {
    ellipse.eccentricity = 0.;
    ellipse.phi = 0.;
    return;
  }
  double mx = 0., my = 0.;
  for (int index : indices) { mx += points[index].x; my += points[index].y; }
  mx /= indices.size(); my /= indices.size();
  double xx = 0., yy = 0., xy = 0.;
  for (int index : indices) {
    double dx = points[index].x - mx, dy = points[index].y - my;
    xx += dx * dx; yy += dy * dy; xy += dx * dy;
  }
  xx /= indices.size(); yy /= indices.size(); xy /= indices.size();
  double phi = 0.5 * std::atan2(2. * xy, xx - yy);
  double c = std::cos(phi), s = std::sin(phi);
  double major = std::max(c * c * xx + 2. * c * s * xy + s * s * yy,
                          s * s * xx - 2. * c * s * xy + c * c * yy);
  double minor = std::min(c * c * xx + 2. * c * s * xy + s * s * yy,
                          s * s * xx - 2. * c * s * xy + c * c * yy);
  if (major <= 0. || minor <= 0.) return;
  if (major < minor) { std::swap(major, minor); phi += 0.5 * M_PI; }
  ellipse.x = mx; ellipse.y = my;
  ellipse.radius = std::sqrt(2. * major);
  ellipse.eccentricity = std::sqrt(std::max(0., 1. - minor / major));
  ellipse.phi = phi;
}

bool has_branch(TTree *tree, const char *name)
{
  if (tree->GetBranch(name))
    return true;
  std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
  return false;
}

bool circle_from_three(const point_t &a, const point_t &b, const point_t &c,
                       circle_t &circle)
{
  const double determinant = 2. * (a.x * (b.y - c.y) +
                                   b.x * (c.y - a.y) +
                                   c.x * (a.y - b.y));
  if (std::abs(determinant) < 1.e-9)
    return false;

  const double aa = a.x * a.x + a.y * a.y;
  const double bb = b.x * b.x + b.y * b.y;
  const double cc = c.x * c.x + c.y * c.y;
  circle.x = (aa * (b.y - c.y) + bb * (c.y - a.y) + cc * (a.y - b.y)) / determinant;
  circle.y = (aa * (c.x - b.x) + bb * (a.x - c.x) + cc * (b.x - a.x)) / determinant;
  circle.radius = std::hypot(a.x - circle.x, a.y - circle.y);
  return std::isfinite(circle.x) && std::isfinite(circle.y) &&
         std::isfinite(circle.radius);
}

double median(std::vector<double> values)
{
  if (values.empty())
    return std::numeric_limits<double>::quiet_NaN();
  auto middle = values.begin() + values.size() / 2;
  std::nth_element(values.begin(), middle, values.end());
  if (values.size() % 2)
    return *middle;
  const double upper = *middle;
  std::nth_element(values.begin(), middle - 1, values.end());
  return 0.5 * (upper + *(middle - 1));
}

circle_t find_one_ring(const std::vector<point_t> &points,
                       const std::vector<int> &indices,
                       std::mt19937 &generator,
                       int iterations,
                       int min_inliers,
                       double tolerance,
                       double time_window,
                       double min_x0, double max_x0,
                       double min_y0, double max_y0,
                       double min_radius, double max_radius,
                       bool force_circle)
{
  circle_t best;
  if (indices.size() < 3)
    return best;

  for (int iteration = 0; iteration < iterations; ++iteration) {
    std::uniform_int_distribution<std::size_t> pick(0, indices.size() - 1);
    const std::size_t ia = pick(generator);

    // Build the minimal temporal hypothesis before fitting the spatial ring.
    // This prevents a triplet from mixing two overlapping rings that are
    // separated in time.
    std::vector<int> temporal_indices;
    const double seed_time = points[indices[ia]].time;
    for (int index : indices)
      if (std::abs(points[index].time - seed_time) <= time_window)
        temporal_indices.push_back(index);
    if (temporal_indices.size() < 3)
      continue;

    std::uniform_int_distribution<std::size_t> temporal_pick(0, temporal_indices.size() - 1);
    const std::size_t ib = temporal_pick(generator);
    const std::size_t ic = temporal_pick(generator);
    if (temporal_indices[ib] == indices[ia] ||
        temporal_indices[ic] == indices[ia] || ib == ic)
      continue;

    circle_t candidate;
    if (!circle_from_three(points[indices[ia]],
                           points[temporal_indices[ib]], points[temporal_indices[ic]], candidate))
      continue;
    if (candidate.x < min_x0 || candidate.x > max_x0 ||
        candidate.y < min_y0 || candidate.y > max_y0 ||
        candidate.radius < min_radius || candidate.radius > max_radius)
      continue;

    std::vector<int> spatial_indices;
    for (int index : indices) {
      const double residual = std::abs(std::hypot(points[index].x - candidate.x,
                                                  points[index].y - candidate.y) -
                                       candidate.radius);
      if (residual <= tolerance) spatial_indices.push_back(index);
    }
    const circle_t circle_candidate = candidate;
    auto collect_inliers = [&](const circle_t &model, double ring_time,
                               std::vector<int> &accepted) {
      accepted.clear();
      for (int index : indices) {
        const double residual = ellipse_residual(model, points[index]);
        if (residual <= tolerance &&
            std::abs(points[index].time - ring_time) <= time_window)
          accepted.push_back(index);
      }
    };

    std::vector<int> accepted_indices;
    std::vector<double> times;
    circle_t circle_model = circle_candidate;
    circle_model.eccentricity = 0.;
    circle_model.phi = 0.;
    circle_model.ring_time = points[indices[ia]].time;
    collect_inliers(circle_model, circle_model.ring_time, accepted_indices);
    times.clear();
    for (int index : accepted_indices)
      times.push_back(points[index].time);
    circle_model.ring_time = median(times);
    collect_inliers(circle_model, circle_model.ring_time, accepted_indices);
    const int circle_inliers = static_cast<int>(accepted_indices.size());
    const std::vector<int> circle_indices = accepted_indices;
    std::vector<double> circle_times;
    for (int index : circle_indices)
      circle_times.push_back(points[index].time);

    candidate = circle_model;
    refine_ellipse(candidate, points, circle_indices, force_circle);
    collect_inliers(candidate, circle_model.ring_time, accepted_indices);
    times.clear();
    for (int index : accepted_indices)
      times.push_back(points[index].time);
    candidate.ring_time = median(times);
    collect_inliers(candidate, candidate.ring_time, accepted_indices);
    const int ellipse_inliers = static_cast<int>(accepted_indices.size());
    if (force_circle ||
        ellipse_inliers < min_inliers || ellipse_inliers <= circle_inliers) {
      candidate = circle_model;
      candidate.ring_time = median(circle_times);
      candidate.inliers = circle_inliers;
      times = circle_times;
    } else {
      candidate.inliers = ellipse_inliers;
    }
    if (static_cast<int>(times.size()) < min_inliers)
      continue;
    int inliers = 0;
    double residual_sum = 0.;
    for (int index : indices) {
      const double spatial = ellipse_residual(candidate, points[index]);
      if (spatial <= tolerance &&
          std::abs(points[index].time - candidate.ring_time) <= time_window) {
        ++inliers;
        residual_sum += spatial;
      }
    }
    if (inliers > best.inliers ||
        (inliers == best.inliers && residual_sum < best.residual_sum)) {
      candidate.inliers = inliers;
      candidate.residual_sum = residual_sum;
      candidate.valid = inliers >= min_inliers;
      best = candidate;
    }
  }
  return best;
}

bool ring_finder(const std::string &filename, const std::string &outfilename,
                 int iterations, int min_inliers, int max_rings,
                 double tolerance, double time_window,
                 double min_x0, double max_x0,
                 double min_y0, double max_y0,
                 double min_radius, double max_radius,
                 bool force_circle)
{
  auto fin = TFile::Open(filename.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << filename << std::endl;
    return false;
  }
  auto frames_in = (TTree *)fin->Get("frames");
  auto cherenkov_in = (TTree *)fin->Get("cherenkov");
  if (!frames_in || !cherenkov_in) {
    std::cerr << "ERROR: input must contain 'frames' and 'cherenkov' trees" << std::endl;
    fin->Close();
    return false;
  }
  for (const char *name : {"nhits", "x", "y", "time"}) {
    if (!has_branch(cherenkov_in, name)) {
      fin->Close();
      return false;
    }
  }
  if (fin->Get("ring")) {
    std::cerr << "ERROR: input already contains a 'ring' tree" << std::endl;
    fin->Close();
    return false;
  }

  const Long64_t entries = frames_in->GetEntries();
  if (cherenkov_in->GetEntries() != entries) {
    std::cerr << "ERROR: frames/cherenkov entry-count mismatch: frames="
              << entries << " cherenkov=" << cherenkov_in->GetEntries() << std::endl;
    fin->Close();
    return false;
  }

  TTreeReader reader(cherenkov_in);
  TTreeReaderValue<UShort_t> ncherenkovhits(reader, "nhits");
  TTreeReaderArray<Float_t> x(reader, "x");
  TTreeReaderArray<Float_t> y(reader, "y");
  TTreeReaderArray<Float_t> time(reader, "time");

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return false;
  }
  fout->cd();
  auto frames_out = frames_in->CloneTree(0);
  auto trigger_in = (TTree *)fin->Get("trigger");
  auto timing_in = (TTree *)fin->Get("timing");
  auto trigger_out = trigger_in ? trigger_in->CloneTree(0) : nullptr;
  auto timing_out = timing_in ? timing_in->CloneTree(0) : nullptr;
  auto cherenkov_out = cherenkov_in->CloneTree(0);
  if (!frames_out || !cherenkov_out || (trigger_in && !trigger_out) ||
      (timing_in && !timing_out)) {
    std::cerr << "ERROR: failed to create output frame trees" << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  UChar_t nring = 0;
  float ring_x0[maxrings];
  float ring_y0[maxrings];
  float ring_r[maxrings];
  float ring_e[maxrings];
  float ring_phi[maxrings];
  float ring_time[maxrings];
  UShort_t ring_ninliers[maxrings];
  auto ring_out = new TTree("ring", "ring candidates, one entry per frame");
  ring_out->Branch("nring", &nring, "nring/b");
  ring_out->Branch("ring_x0", ring_x0, "ring_x0[nring]/F");
  ring_out->Branch("ring_y0", ring_y0, "ring_y0[nring]/F");
  ring_out->Branch("ring_r", ring_r, "ring_r[nring]/F");
  ring_out->Branch("ring_e", ring_e, "ring_e[nring]/F");
  ring_out->Branch("ring_phi", ring_phi, "ring_phi[nring]/F");
  ring_out->Branch("ring_time", ring_time, "ring_time[nring]/F");
  ring_out->Branch("ring_ninliers", ring_ninliers, "ring_ninliers[nring]/s");

  std::mt19937 generator(0x9e3779b9u);
  Long64_t frames = 0;
  Long64_t accepted = 0;
  for (Long64_t iframe = 0; iframe < entries; ++iframe) {
    if (frames_in->GetEntry(iframe) <= 0 || !reader.Next()) {
      std::cerr << "ERROR: failed to read synchronized frame entry " << iframe << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    ++frames;
    nring = 0;
    const int nhits = static_cast<int>(*ncherenkovhits);
    if (x.GetSize() < nhits || y.GetSize() < nhits || time.GetSize() < nhits) {
      std::cerr << "ERROR: invalid Cherenkov hit array at frame entry " << iframe
                << " nhits=" << nhits << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    std::vector<point_t> points;
    std::vector<int> remaining;
    for (int index = 0; index < nhits; ++index) {
      if (std::isfinite(x[index]) && std::isfinite(y[index]) &&
          std::isfinite(time[index])) {
        points.push_back({x[index], y[index], time[index]});
        remaining.push_back(static_cast<int>(points.size() - 1));
      }
    }

    while (nring < max_rings && remaining.size() >= 3) {
      circle_t ring = find_one_ring(points, remaining, generator, iterations,
                                    min_inliers, tolerance, time_window,
                                    min_x0, max_x0, min_y0, max_y0,
                                    min_radius, max_radius, force_circle);
      if (!ring.valid)
        break;
      ring_x0[nring] = static_cast<float>(ring.x);
      ring_y0[nring] = static_cast<float>(ring.y);
      ring_r[nring] = static_cast<float>(ring.radius);
      ring_e[nring] = static_cast<float>(ring.eccentricity);
      ring_phi[nring] = static_cast<float>(ring.phi);
      ring_time[nring] = static_cast<float>(ring.ring_time);
      ring_ninliers[nring] = static_cast<UShort_t>(ring.inliers);
      ++nring;
      ++accepted;

      std::vector<int> next;
      for (int index : remaining) {
        const double spatial = ellipse_residual(ring, points[index]);
        if (spatial > tolerance ||
            std::abs(points[index].time - ring.ring_time) > time_window)
          next.push_back(index);
      }
      if (next.size() == remaining.size())
        break;
      remaining.swap(next);
    }

    if (trigger_in)
      trigger_in->GetEntry(iframe);
    if (timing_in)
      timing_in->GetEntry(iframe);
    cherenkov_in->GetEntry(iframe);
    frames_out->Fill();
    if (trigger_out)
      trigger_out->Fill();
    if (timing_out)
      timing_out->Fill();
    cherenkov_out->Fill();
    ring_out->Fill();
  }

  for (auto tree : {frames_out, trigger_out, timing_out, cherenkov_out, ring_out}) {
    if (tree && tree->GetEntries() != entries) {
      std::cerr << "ERROR: output tree entry-count mismatch: " << tree->GetName()
                << "=" << tree->GetEntries() << " expected=" << entries << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }
  }
  auto participation = (TTree *)fin->Get("spill_participation");
  if (participation) {
    fout->cd();
    if (!participation->CloneTree(-1, "fast")) {
      std::cerr << "ERROR: failed to copy spill_participation tree" << std::endl;
      return false;
    }
  }
  fout->Write();
  fout->Close();
  fin->Close();
  std::cout << "frames processed: " << frames << std::endl
            << "rings found:      " << accepted << std::endl
            << "output:            " << outfilename << std::endl;
  return true;
}

}

int main(int argc, char **argv)
{
  namespace po = boost::program_options;
  std::string input, output;
  int iterations = 512, min_inliers = 8, max_rings = maxrings;
  double tolerance = 3.5, time_window = 5.;
  double min_x0 = -100., max_x0 = 100.;
  double min_y0 = -100., max_y0 = 100.;
  double min_radius = 1., max_radius = 200.;
  bool force_circle = false;
  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input triggered ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("iterations", po::value<int>(&iterations)->default_value(iterations), "RANSAC iterations per ring")
    ("min-inliers", po::value<int>(&min_inliers)->default_value(min_inliers), "minimum spatial/time inliers")
    ("max-rings", po::value<int>(&max_rings)->default_value(max_rings), "maximum rings per frame")
    ("tolerance", po::value<double>(&tolerance)->default_value(tolerance), "radial tolerance")
    ("time-window", po::value<double>(&time_window)->default_value(time_window), "ring hit time window")
    ("min-x0", po::value<double>(&min_x0)->default_value(min_x0), "minimum ring center x")
    ("max-x0", po::value<double>(&max_x0)->default_value(max_x0), "maximum ring center x")
    ("min-y0", po::value<double>(&min_y0)->default_value(min_y0), "minimum ring center y")
    ("max-y0", po::value<double>(&max_y0)->default_value(max_y0), "maximum ring center y")
    ("min-radius", po::value<double>(&min_radius)->default_value(min_radius), "minimum ring radius")
    ("max-radius", po::value<double>(&max_radius)->default_value(max_radius), "maximum ring radius")
    ("force-circle", po::bool_switch(&force_circle), "force eccentricity to zero")
    ;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) { std::cout << options << std::endl; return 0; }
    po::notify(vm);
    if (iterations < 1 || min_inliers < 3 || max_rings < 1 ||
        max_rings > maxrings ||
        tolerance <= 0. || time_window < 0. || min_x0 > max_x0 ||
        min_y0 > max_y0 || min_radius <= 0. || min_radius > max_radius)
      throw std::runtime_error("invalid ring-finder parameter");
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << std::endl << options << std::endl;
    return 1;
  }
  return ring_finder(input, output, iterations, min_inliers, max_rings,
                     tolerance, time_window, min_x0, max_x0, min_y0, max_y0,
                     min_radius, max_radius, force_circle) ? 0 : 1;
}

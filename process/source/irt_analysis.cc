#include "trigger_reader.h"
#include <boost/program_options.hpp>
#include <TFile.h>
#include <TH2D.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <string>
#ifdef IRT_ANALYSIS_HAS_CUDA
#include "irt_analysis_cuda.h"
#endif

namespace po = boost::program_options;

int main(int argc, char **argv)
{
  std::string input, output, ring_name = "ring";
  long long max_events = -1;
  double x0min = -INFINITY, x0max = INFINITY, y0min = -INFINITY, y0max = INFINITY;
  double rmin = 0., rmax = INFINITY;
  long long min_hits = 0, max_hits = std::numeric_limits<long long>::max();
  po::options_description options("options");
  options.add_options()
    ("help,h", "show help")
    ("input", po::value<std::string>(&input)->required(), "IRT ROOT input")
    ("output", po::value<std::string>(&output)->required(), "ROOT output")
    ("ring", po::value<std::string>(&ring_name)->default_value(ring_name), "ring tree name")
    ("max-events", po::value<long long>(&max_events)->default_value(max_events), "maximum frames")
    ("ring-x0-min", po::value<double>(&x0min)->default_value(x0min), "minimum ring x0")
    ("ring-x0-max", po::value<double>(&x0max)->default_value(x0max), "maximum ring x0")
    ("ring-y0-min", po::value<double>(&y0min)->default_value(y0min), "minimum ring y0")
    ("ring-y0-max", po::value<double>(&y0max)->default_value(y0max), "maximum ring y0")
    ("ring-r-min", po::value<double>(&rmin)->default_value(rmin), "minimum ring radius")
    ("ring-r-max", po::value<double>(&rmax)->default_value(rmax), "maximum ring radius")
    ("min-cherenkov-hits", po::value<long long>(&min_hits)->default_value(min_hits), "minimum Cherenkov hits")
    ("max-cherenkov-hits", po::value<long long>(&max_hits)->default_value(max_hits), "maximum Cherenkov hits")
    ("gpu", po::bool_switch(), "use CUDA for the theta-time accumulator");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) { std::cout << options << '\n'; return 0; }
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << '\n' << options << '\n'; return 1;
  }
  trigger_reader_t reader;
  if (!reader.open(input, ring_name)) return 1;
  if (!reader.has_rings()) { std::cerr << "ERROR: missing ring tree\n"; return 1; }
#ifdef IRT_ANALYSIS_HAS_CUDA
  const bool use_gpu = vm["gpu"].as<bool>();
  if (use_gpu && !irt_analysis_cuda::available()) { std::cerr << "ERROR: CUDA requested but unavailable\n"; return 1; }
#else
  const bool use_gpu = false;
  if (vm["gpu"].as<bool>()) { std::cerr << "ERROR: this build has no CUDA backend\n"; return 1; }
#endif
  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) { std::cerr << "ERROR: cannot create output\n"; return 1; }
  constexpr int nt = 100, ntime = 64;
  constexpr double tmin = 0., tmax = .1, timemin = -32., timemax = 32.;
  constexpr double rt = .002, rtime = 1., cut = 4.;
  TH2D hist("hAverageCherenkovThetaVsNHits", "average Cherenkov angle;selected hits;#theta [rad]", 128, 0., 128., 512, 0., .1);
  long long frames = 0, one_ring = 0, filled = 0;
  while (reader.next_spill() && (max_events < 0 || frames < max_events)) {
    while (reader.next_frame() && (max_events < 0 || frames < max_events)) {
      ++frames;
      if (reader.rings().size() != 1) continue;
      const auto &ring = reader.rings().front();
      if (ring.x0 < x0min || ring.x0 > x0max || ring.y0 < y0min ||
          ring.y0 > y0max || ring.radius < rmin || ring.radius > rmax)
        continue;
      const auto nhits = static_cast<long long>(reader.cherenkov_hits().size());
      if (nhits < min_hits || nhits > max_hits)
        continue;
      ++one_ring;
      std::vector<const hit_t *> hits;
      for (const auto &h : reader.cherenkov_hits())
        if (std::isfinite(h.theta) && std::isfinite(h.time) && h.theta >= tmin && h.theta <= tmax && h.time >= timemin && h.time <= timemax)
          hits.push_back(&h);
      double score_max = -1., best_theta = 0., best_time = 0.;
#ifdef IRT_ANALYSIS_HAS_CUDA
      if (use_gpu) {
        std::vector<float> theta_gpu, time_gpu;
        theta_gpu.reserve(hits.size()); time_gpu.reserve(hits.size());
        for (const auto *h : hits) { theta_gpu.push_back(h->theta); time_gpu.push_back(h->time); }
        irt_analysis_cuda::peak_t peak;
        if (!irt_analysis_cuda::peak(theta_gpu, time_gpu, peak)) continue;
        best_theta = peak.theta; best_time = peak.time;
      } else
#endif
      for (int it = 0; it < nt; ++it) for (int iu = 0; iu < ntime; ++iu) {
        const double theta = tmin + (it + .5) * (tmax-tmin)/nt;
        const double time = timemin + (iu + .5) * (timemax-timemin)/ntime;
        double score = 0.;
        for (const auto *h : hits) { const double a=(h->theta-theta)/rt, b=(h->time-time)/rtime; score += std::exp(-.5*(a*a+b*b)); }
        if (score > score_max) { score_max=score; best_theta=theta; best_time=time; }
      }
      double sum=0.; int n=0;
      for (const auto *h : hits)
        if (std::abs(h->theta-best_theta) <= cut*rt && std::abs(h->time-best_time) <= cut*rtime) { sum += h->theta; ++n; }
      if (n > 0) { hist.Fill(n, sum/n); ++filled; }
    }
  }
  file->cd(); hist.Write(); file->Close();
  std::cout << "frames processed: " << frames << '\n'
            << "one-ring frames: " << one_ring << '\n'
            << "histogram fills: " << filled << '\n'
            << "output: " << output << '\n';
  return 0;
}

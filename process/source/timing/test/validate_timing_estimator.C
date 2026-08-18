#include "../src/TimingEstimator.cxx"

#include <TSystem.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

struct diff_acc_t {
  double max_abs = 0.;
  double sum2 = 0.;

  void add(double value)
  {
    double a = std::abs(value);
    if (a > max_abs)
      max_abs = a;
    sum2 += value * value;
  }

  double rms(int n) const
  {
    return n > 0 ? std::sqrt(sum2 / n) : 0.;
  }
};

void validate_timing_estimator(int nevents = 1000,
                               const char *model_dir = "../python",
                               const char *python = "python3")
{
  std::string macro_dir = gSystem->DirName(__FILE__);
  std::string model_path = model_dir;
  if (!model_path.empty() && model_path[0] != '/')
    model_path = macro_dir + "/" + model_path;

  std::string prefix = std::string("/tmp/timing_estimator_validate_") +
                       std::to_string(gSystem->GetPid());
  std::string events_path = prefix + "_events.csv";
  std::string ref_path = prefix + "_python.csv";

  std::vector<std::array<double, 64>> events;
  events.reserve(nevents);

  std::mt19937_64 rng(0x20260818);
  std::uniform_real_distribution<double> event_time(-200., 200.);
  std::normal_distribution<double> jitter(0., 0.35);
  std::normal_distribution<double> shape(0., 0.8);

  std::ofstream events_file(events_path);
  events_file << std::setprecision(17);
  for (int iev = 0; iev < nevents; ++iev) {
    std::array<double, 64> event;
    double t0 = event_time(rng);
    double t1 = t0 + 0.1 * jitter(rng);
    for (int i = 0; i < 32; ++i)
      event[i] = t0 + 0.03 * i + jitter(rng) + 0.2 * shape(rng);
    for (int i = 0; i < 32; ++i)
      event[32 + i] = t1 - 0.02 * i + jitter(rng) + 0.2 * shape(rng);

    for (int i = 0; i < 64; ++i) {
      if (i)
        events_file << ",";
      events_file << event[i];
    }
    events_file << "\n";
    events.push_back(event);
  }
  events_file.close();

  std::string helper = macro_dir + "/python_reference.py";
  std::ostringstream cmd;
  if (const char *ldpath = std::getenv("TIMING_ESTIMATOR_LD_LIBRARY_PATH")) {
    cmd << "LD_LIBRARY_PATH=" << ldpath << " ";
  }
  if (const char *pythonpath = std::getenv("TIMING_ESTIMATOR_PYTHONPATH")) {
    cmd << "PYTHONPATH=" << pythonpath << " ";
  }
  cmd << python
      << " " << helper
      << " --model-dir " << model_path
      << " --input " << events_path
      << " --output " << ref_path;

  int ret = gSystem->Exec(cmd.str().c_str());
  if (ret != 0) {
    std::cerr << "ERROR: Python reference failed with status " << ret << std::endl;
    return;
  }

  std::ifstream ref(ref_path);
  std::string line;
  std::getline(ref, line);

  TimingEstimator estimator;
  diff_acc_t dT0, ds0, dT1, ds1, dT, dsT, dDeltaT;
  int nread = 0;

  while (std::getline(ref, line)) {
    if (line.empty())
      continue;

    std::stringstream ss(line);
    std::string item;
    std::array<double, 7> py;
    for (int i = 0; i < 7; ++i) {
      if (!std::getline(ss, item, ',')) {
        std::cerr << "ERROR: malformed Python reference row " << nread << std::endl;
        return;
      }
      py[i] = std::stod(item);
    }

    auto cxx = estimator.estimate(events[nread]);
    double cxx_delta = cxx.T1 - cxx.T0;

    dT0.add(cxx.T0 - py[0]);
    ds0.add(cxx.sigma0 - py[1]);
    dT1.add(cxx.T1 - py[2]);
    ds1.add(cxx.sigma1 - py[3]);
    dT.add(cxx.T - py[4]);
    dsT.add(cxx.sigmaT - py[5]);
    dDeltaT.add(cxx_delta - py[6]);
    ++nread;
  }

  auto print = [&](const char *name, const diff_acc_t &diff) {
    std::cout << std::setw(14) << name
              << " max_abs=" << std::setprecision(12) << diff.max_abs
              << " rms=" << diff.rms(nread)
              << std::endl;
  };

  std::cout << "tested events: " << nread << std::endl;
  print("T0", dT0);
  print("sigma0", ds0);
  print("T1", dT1);
  print("sigma1", ds1);
  print("T", dT);
  print("sigmaT", dsT);
  print("DeltaT", dDeltaT);

  std::array<double, 64> shifted = events.front();
  auto ref0 = estimator.estimate(shifted);
  constexpr double delta0 = 17.25;
  constexpr double delta1 = -9.75;
  for (int i = 0; i < 32; ++i)
    shifted[i] += delta0;
  for (int i = 32; i < 64; ++i)
    shifted[i] += delta1;
  auto sh = estimator.estimate(shifted);

  std::cout << "translation check:" << std::endl;
  std::cout << "  T0 shift      = " << sh.T0 - ref0.T0
            << " expected=" << delta0 << std::endl;
  std::cout << "  T1 shift      = " << sh.T1 - ref0.T1
            << " expected=" << delta1 << std::endl;
  std::cout << "  sigma0 change = " << sh.sigma0 - ref0.sigma0 << std::endl;
  std::cout << "  sigma1 change = " << sh.sigma1 - ref0.sigma1 << std::endl;

  gSystem->Unlink(events_path.c_str());
  gSystem->Unlink(ref_path.c_str());
}

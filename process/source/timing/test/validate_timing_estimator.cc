#include "../include/TimingEstimator.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

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

int
main(int argc, char **argv)
{
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " events.csv python_reference.csv" << std::endl;
    return 1;
  }

  std::ifstream events_file(argv[1]);
  std::ifstream ref_file(argv[2]);
  if (!events_file) {
    std::cerr << "ERROR: could not open events file: " << argv[1] << std::endl;
    return 1;
  }
  if (!ref_file) {
    std::cerr << "ERROR: could not open Python reference file: " << argv[2] << std::endl;
    return 1;
  }

  std::string line;
  std::getline(ref_file, line);

  TimingEstimator estimator;
  diff_acc_t dT0, ds0, dT1, ds1, dT, dsT, dDeltaT;
  int nread = 0;

  while (std::getline(events_file, line)) {
    if (line.empty())
      continue;

    std::array<double, 64> event;
    std::stringstream es(line);
    std::string item;
    for (int i = 0; i < 64; ++i) {
      if (!std::getline(es, item, ',')) {
        std::cerr << "ERROR: malformed event row " << nread << std::endl;
        return 1;
      }
      event[i] = std::stod(item);
    }

    if (!std::getline(ref_file, line)) {
      std::cerr << "ERROR: missing Python reference row " << nread << std::endl;
      return 1;
    }

    std::array<double, 7> py;
    std::stringstream rs(line);
    for (int i = 0; i < 7; ++i) {
      if (!std::getline(rs, item, ',')) {
        std::cerr << "ERROR: malformed Python reference row " << nread << std::endl;
        return 1;
      }
      py[i] = std::stod(item);
    }

    auto cxx = estimator.estimate(event);
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

  std::array<double, 64> shifted;
  for (int i = 0; i < 64; ++i)
    shifted[i] = 100. + 0.1 * i;
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

  return 0;
}

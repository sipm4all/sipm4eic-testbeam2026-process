#include <boost/program_options.hpp>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <limits>

namespace po = boost::program_options;

int main(int argc, char **argv)
{
  std::string input, output;
  int events = 1000, hits = 32;
  unsigned int seed = 12345;
  double x0_override = std::numeric_limits<double>::quiet_NaN();
  double y0_override = std::numeric_limits<double>::quiet_NaN();
  double a_override = std::numeric_limits<double>::quiet_NaN();
  double b_override = std::numeric_limits<double>::quiet_NaN();
  double theta_override = std::numeric_limits<double>::quiet_NaN();
  double sigma_override = std::numeric_limits<double>::quiet_NaN();
  po::options_description options("options");
  options.add_options()
    ("help,h", "show help")
    ("input", po::value<std::string>(&input)->required(), "ROOT file containing hResults")
    ("output", po::value<std::string>(&output)->required(), "synthetic triggered ROOT file")
    ("events", po::value<int>(&events)->default_value(events), "number of generated frames")
    ("hits", po::value<int>(&hits)->default_value(hits), "Cherenkov hits per frame")
    ("seed", po::value<unsigned int>(&seed)->default_value(seed), "random seed")
    ("x0", po::value<double>(&x0_override), "override ellipse X0")
    ("y0", po::value<double>(&y0_override), "override ellipse Y0")
    ("a", po::value<double>(&a_override), "override ellipse semi-axis A")
    ("b", po::value<double>(&b_override), "override ellipse semi-axis B")
    ("theta", po::value<double>(&theta_override), "override ellipse angle")
    ("sigma-rho", po::value<double>(&sigma_override), "override normalized radial width");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) { std::cout << options << '\n'; return 0; }
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << '\n' << options << '\n'; return 1;
  }
  if (events <= 0 || hits <= 0 || hits > 65535) {
    std::cerr << "ERROR: invalid events or hits\n"; return 1;
  }

  TFile source(input.c_str(), "READ");
  auto *results = dynamic_cast<TH1 *>(source.Get("hResults"));
  if (!results) { std::cerr << "ERROR: missing hResults in " << input << '\n'; return 1; }
  auto value = [&](const char *label) {
    return results->GetBinContent(results->GetXaxis()->FindBin(label));
  };
  const double x0 = std::isnan(x0_override) ? value("X0") : x0_override;
  const double y0 = std::isnan(y0_override) ? value("Y0") : y0_override;
  const double a = std::isnan(a_override) ? value("A") : a_override;
  const double b = std::isnan(b_override) ? value("B") : b_override;
  const double angle = std::isnan(theta_override) ? value("theta") : theta_override;
  const double sigma_rho = std::isnan(sigma_override) ? value("sigmaRho") : sigma_override;
  if (!(a > 0. && b > 0.) || !std::isfinite(sigma_rho) || sigma_rho < 0.) {
    std::cerr << "ERROR: invalid ellipse parameters in hResults\n"; return 1;
  }

  TFile destination(output.c_str(), "RECREATE");
  if (destination.IsZombie()) { std::cerr << "ERROR: could not create " << output << '\n'; return 1; }
  Int_t spill = 0, frame = 0;
  TTree frames("frames", "synthetic frame index");
  frames.Branch("spill", &spill, "spill/I");
  frames.Branch("frame", &frame, "frame/I");
  UShort_t nhits = 0;
  Float_t x[65535] = {}, y[65535] = {}, time[65535] = {};
  TTree cherenkov("cherenkov", "synthetic Cherenkov hits");
  cherenkov.Branch("nhits", &nhits, "nhits/s");
  cherenkov.Branch("x", x, "x[nhits]/F");
  cherenkov.Branch("y", y, "y[nhits]/F");
  cherenkov.Branch("time", time, "time[nhits]/F");
  UChar_t nrings = 0;
  Float_t ring_x0[1] = {}, ring_y0[1] = {}, ring_r[1] = {}, ring_e[1] = {};
  Float_t ring_phi[1] = {}, ring_time[1] = {};
  UShort_t ring_ninliers[1] = {};
  TTree ring("ring", "synthetic ellipse candidates");
  ring.Branch("nrings", &nrings, "nrings/b");
  ring.Branch("x0", ring_x0, "x0[nrings]/F");
  ring.Branch("y0", ring_y0, "y0[nrings]/F");
  ring.Branch("r", ring_r, "r[nrings]/F");
  ring.Branch("e", ring_e, "e[nrings]/F");
  ring.Branch("phi", ring_phi, "phi[nrings]/F");
  ring.Branch("time", ring_time, "time[nrings]/F");
  ring.Branch("ninliers", ring_ninliers, "ninliers[nrings]/s");

  std::mt19937 generator(seed);
  std::normal_distribution<double> smear(0., sigma_rho);
  std::uniform_real_distribution<double> azimuth(0., 2. * 3.14159265358979323846);
  for (int iev = 0; iev < events; ++iev) {
    frame = iev; nhits = static_cast<UShort_t>(hits);
    nrings = 1; ring_x0[0] = x0; ring_y0[0] = y0;
    ring_r[0] = static_cast<Float_t>(std::max(a, b));
    ring_e[0] = static_cast<Float_t>(std::sqrt(1. -
      (std::min(a, b) * std::min(a, b)) / (std::max(a, b) * std::max(a, b))));
    ring_phi[0] = static_cast<Float_t>(angle); ring_time[0] = 0.f;
    ring_ninliers[0] = static_cast<UShort_t>(hits);
    for (int ihit = 0; ihit < hits; ++ihit) {
      const double t = azimuth(generator);
      const double q = smear(generator);
      const double u0 = a * std::cos(t), v0 = b * std::sin(t);
      const double u = u0 * (1. + q);
      const double v = v0 * (1. + q);
      x[ihit] = static_cast<Float_t>(x0 + u * std::cos(angle) - v * std::sin(angle));
      y[ihit] = static_cast<Float_t>(y0 + u * std::sin(angle) + v * std::cos(angle));
      time[ihit] = 0.f;
    }
    frames.Fill(); cherenkov.Fill(); ring.Fill();
  }
  frames.Write(); cherenkov.Write(); ring.Write(); destination.Close();
  std::cout << "ellipse parameters: X0=" << x0 << " Y0=" << y0
            << " A=" << a << " B=" << b << " theta=" << angle
            << " sigmaRho=" << sigma_rho << '\n'
            << "generated events: " << events << '\n'
            << "generated hits: " << static_cast<long long>(events) * hits << '\n'
            << "output: " << output << '\n';
  return 0;
}

#include <TFile.h>
#include <TH1.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct ellipse_fit_t {
  double nsig = 0.;
  double x0 = 0.;
  double y0 = 0.;
  double a = 0.;
  double b = 0.;
  double angle = 0.;
  double sigma_rho = 0.;
};

double parameter(const TH1 &results, const char *name)
{
  for (int bin = 1; bin <= results.GetNbinsX(); ++bin) {
    if (std::string(results.GetXaxis()->GetBinLabel(bin)) == name)
      return results.GetBinContent(bin);
  }
  throw std::runtime_error(std::string("missing fit parameter: ") + name);
}

double signal_density(const ellipse_fit_t &fit, double x, double y,
                      double n_sigma)
{
  const double dx = x - fit.x0;
  const double dy = y - fit.y0;
  const double cosine = std::cos(fit.angle);
  const double sine = std::sin(fit.angle);
  const double u = (dx * cosine + dy * sine) / fit.a;
  const double v = (-dx * sine + dy * cosine) / fit.b;
  const double rho = std::hypot(u, v);
  if (std::abs(rho - 1.) > n_sigma * fit.sigma_rho)
    return 0.;

  return fit.nsig / (2. * TMath::Pi() * fit.a * fit.b) *
         TMath::Gaus(rho, 1., fit.sigma_rho, true);
}

double integrate_pixels(const ellipse_fit_t &fit,
                        const std::vector<std::pair<double, double>> &centres,
                        double pixel_size, double integration_step,
                        double n_sigma)
{
  const int samples = std::max(
      1, static_cast<int>(std::ceil(pixel_size / integration_step)));
  const double step = pixel_size / samples;
  double integral = 0.;
  for (const auto &[center_x, center_y] : centres) {
    for (int ix = 0; ix < samples; ++ix) {
      const double x = center_x - pixel_size / 2. + (ix + 0.5) * step;
      for (int iy = 0; iy < samples; ++iy) {
        const double y = center_y - pixel_size / 2. + (iy + 0.5) * step;
        integral += signal_density(fit, x, y, n_sigma) * step * step;
      }
    }
  }
  return integral;
}

std::vector<std::pair<double, double>> ideal_lattice(
    const ellipse_fit_t &fit, double pixel_size, double pixel_pitch,
    double n_sigma, double phase_x, double phase_y)
{
  const double extent = std::max(std::abs(fit.x0), std::abs(fit.y0)) +
      std::max(fit.a, fit.b) * (1. + n_sigma * fit.sigma_rho) +
      pixel_size + pixel_pitch;
  std::vector<std::pair<double, double>> centres;
  for (double x = -extent + phase_x; x <= extent; x += pixel_pitch)
    for (double y = -extent + phase_y; y <= extent; y += pixel_pitch)
      centres.emplace_back(x, y);
  return centres;
}

std::vector<std::pair<double, double>> complete_pdu_geometry(
    double pixel_size, double pixel_pitch, double pdu_gap)
{
  constexpr int pixels_per_side = 16;
  const double pdu_extent = pixel_size + (pixels_per_side - 1) * pixel_pitch;
  const double pdu_pitch = pdu_extent + pdu_gap;
  std::vector<std::pair<double, double>> centres;
  for (int pdu_x = -1; pdu_x <= 1; ++pdu_x) {
    for (int pdu_y = -1; pdu_y <= 1; ++pdu_y) {
      for (int column = 0; column < pixels_per_side; ++column) {
        for (int row = 0; row < pixels_per_side; ++row) {
          centres.emplace_back(
              pdu_x * pdu_pitch + (column - 7.5) * pixel_pitch,
              pdu_y * pdu_pitch + (row - 7.5) * pixel_pitch);
        }
      }
    }
  }
  return centres;
}

std::vector<std::pair<double, double>> real_geometry(double pixel_pitch)
{
  // Bottom-left PDU placements used by process/source/geometry.h.
  const std::vector<std::pair<double, double>> placements = {
      {-82., 30.}, {-26., 35.}, {30., 30.},
      {-82., -26.},             {30., -26.},
      {-82., -82.}, {-26., -87.}, {30., -82.}};
  std::vector<std::pair<double, double>> centres;
  for (const auto &[placement_x, placement_y] : placements) {
    for (int column = 0; column < 16; ++column) {
      for (int row = 0; row < 16; ++row) {
        // The 1.85 mm offset and 0.3 mm matrix gaps reproduce geometry.h.
        const double x = placement_x + 1.85 + pixel_pitch * column +
                         (column > 7 ? 0.3 : 0.);
        const double y = placement_y + 1.85 + pixel_pitch * row +
                         (row > 7 ? 0.3 : 0.);
        centres.emplace_back(x, y);
      }
    }
  }
  return centres;
}

} // namespace

void hitmap_acceptance(const char *fit_filename,
                       double n_sigma = 3.,
                       double pixel_size = 3.,
                       double pixel_pitch = 3.2,
                       double pdu_gap = 3.,
                       double integration_step = 0.05,
                       int lattice_phase_steps = 8)
{
  TFile input(fit_filename, "READ");
  if (input.IsZombie()) {
    std::cerr << "ERROR: could not open fit file: " << fit_filename << '\n';
    return;
  }
  const auto *results = dynamic_cast<TH1 *>(input.Get("hResults"));
  if (!results) {
    std::cerr << "ERROR: input does not contain hResults\n";
    return;
  }

  ellipse_fit_t fit;
  try {
    fit = {parameter(*results, "Nsig"), parameter(*results, "X0"),
           parameter(*results, "Y0"), parameter(*results, "A"),
           parameter(*results, "B"), parameter(*results, "theta"),
           parameter(*results, "sigmaRho")};
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return;
  }
  if (fit.nsig <= 0. || fit.a <= 0. || fit.b <= 0. ||
      fit.sigma_rho <= 0. || n_sigma <= 0. || pixel_size <= 0. ||
      pixel_pitch <= 0. || integration_step <= 0. ||
      lattice_phase_steps <= 0) {
    std::cerr << "ERROR: invalid fit or integration parameters\n";
    return;
  }

  // The rho Jacobian is proportional to rho. In a symmetric interval around
  // rho=1 its odd contribution cancels, leaving the Gaussian containment.
  const double continuous = fit.nsig * TMath::Erf(n_sigma / std::sqrt(2.));

  double ideal_sum = 0.;
  double ideal_min = std::numeric_limits<double>::max();
  double ideal_max = 0.;
  for (int ix = 0; ix < lattice_phase_steps; ++ix) {
    for (int iy = 0; iy < lattice_phase_steps; ++iy) {
      const double ideal = integrate_pixels(
          fit,
          ideal_lattice(fit, pixel_size, pixel_pitch, n_sigma,
                        ix * pixel_pitch / lattice_phase_steps,
                        iy * pixel_pitch / lattice_phase_steps),
          pixel_size, integration_step, n_sigma);
      ideal_sum += ideal;
      ideal_min = std::min(ideal_min, ideal);
      ideal_max = std::max(ideal_max, ideal);
    }
  }
  const double ideal = ideal_sum /
      (lattice_phase_steps * lattice_phase_steps);
  const double complete_pdu = integrate_pixels(
      fit, complete_pdu_geometry(pixel_size, pixel_pitch, pdu_gap),
      pixel_size, integration_step, n_sigma);
  const double real = integrate_pixels(
      fit, real_geometry(pixel_pitch), pixel_size, integration_step, n_sigma);

  auto percentage = [](double numerator, double denominator) {
    return denominator > 0. ? 100. * numerator / denominator : 0.;
  };
  std::cout << "fit: " << fit_filename
            << "\nring interval: +/- " << n_sigma << " sigmaRho"
            << "\npixel: " << pixel_size << " x " << pixel_size
            << " mm2, pitch: " << pixel_pitch << " mm"
            << "\nPDU gap: " << pdu_gap << " mm"
            << "\ncontinuous ring integral: " << continuous
            << "\nideal full-array integral: " << ideal
            << " (phase range " << ideal_min << " to " << ideal_max << ")"
            << "\ncomplete 3x3-PDU integral: " << complete_pdu
            << "\nreal eight-PDU integral: " << real
            << "\nideal/continuous: " << percentage(ideal, continuous) << " %"
            << "\ncomplete-PDU/continuous: "
            << percentage(complete_pdu, continuous) << " %"
            << "\ncomplete-PDU/ideal: " << percentage(complete_pdu, ideal) << " %"
            << "\nreal/continuous: " << percentage(real, continuous) << " %"
            << "\nreal/ideal: " << percentage(real, ideal) << " %"
            << "\nreal/complete-PDU: " << percentage(real, complete_pdu) << " %\n";
}

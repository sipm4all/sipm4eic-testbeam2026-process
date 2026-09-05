#include "irt_geometry.h"

#include "trigger_reader.h"

#include <boost/program_options.hpp>
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace po = boost::program_options;

namespace {

bool read_geometry_config(const std::string &filename, irt::geometry_t &geometry)
{
  std::ifstream input(filename);
  if (!input) return false;
  std::string section, line;
  auto values = [](const std::string &text, double out[3]) {
    std::string cleaned = text;
    for (char &c : cleaned)
      if (c == '(' || c == ')' || c == ',') c = ' ';
    std::istringstream stream(cleaned);
    return static_cast<bool>(stream >> out[0] >> out[1] >> out[2]);
  };
  while (std::getline(input, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream stream(line);
    std::string key, equals;
    if (!(stream >> key)) continue;
    if (key.front() == '[') { section = key; continue; }
    stream >> equals;
    std::string rest;
    std::getline(stream, rest);
    double v[3] = {};
    if (section == "[mirror]" && key == "center" && values(rest, v))
      geometry.mirror_center.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[mirror]" && key == "radius")
      geometry.mirror_radius = std::stod(rest);
    else if (section == "[mirror]" && key == "pivot_center" && values(rest, v))
      geometry.mirror_pivot.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[mirror]" && key == "rotation_vector" && values(rest, v))
      geometry.mirror_rotation_vector.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[detector]" && key == "sphere_center" && values(rest, v))
      geometry.detector_center.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[detector]" && key == "sphere_radius")
      geometry.detector_radius = std::stod(rest);
    else if (section == "[detector]" && key == "theta")
      geometry.detector_theta = std::stod(rest);
    else if (section == "[detector]" && key == "phi")
      geometry.detector_phi = std::stod(rest);
    else if (section == "[detector]" && key == "rotation_vector" && values(rest, v))
      geometry.detector_rotation_vector.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[assembly]" && key == "rotation_x_pivot" && values(rest, v))
      geometry.assembly_rotation_x_pivot.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[assembly]" && key == "rotation_y_pivot" && values(rest, v))
      geometry.assembly_rotation_y_pivot.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[assembly]" && key == "rotation_x")
      geometry.assembly_rotation_x = std::stod(rest);
    else if (section == "[assembly]" && key == "rotation_y")
      geometry.assembly_rotation_y = std::stod(rest);
    else if (section == "[track]" && key == "origin" && values(rest, v))
      geometry.track_origin.SetXYZ(v[0], v[1], v[2]);
    else if (section == "[track]" && key == "eta")
      geometry.track_eta = std::stod(rest);
    else if (section == "[track]" && key == "phi")
      geometry.track_phi = std::stod(rest);
    else if (section == "[track]" && key == "emission_z")
      geometry.emission_z = std::stod(rest);
  }
  return geometry.mirror_radius > 0. && geometry.detector_radius > 0.;
}

struct fit_parameter_t { std::string name; double start, min, max; bool free; };
bool read_fit_config(const std::string &filename, std::vector<fit_parameter_t> &out)
{
  std::ifstream input(filename);
  if (!input) return false;
  std::string line;
  while (std::getline(input, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream stream(line);
    fit_parameter_t p; std::string status;
    if (!(stream >> p.name)) continue;
    if (!(stream >> p.start >> p.min >> p.max >> status) ||
        (status != "free" && status != "fixed") || p.min > p.max)
      return false;
    p.free = status == "free";
    out.push_back(p);
  }
  return !out.empty();
}

void apply_fit_parameters(irt::geometry_t &g,
                          const std::vector<fit_parameter_t> &parameters,
                          const double *values)
{
  TVector3 mirror_delta, mirror_rotation, detector_delta, detector_rotation;
  double eta = g.track_eta, phi = g.track_phi, emission_z = g.emission_z;
  TVector3 origin = g.track_origin;
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    const auto &n = parameters[i].name; const double v = values[i];
    if (n == "track_eta") eta = v; else if (n == "track_phi") phi = v;
    else if (n == "track_origin_x") origin.SetX(v);
    else if (n == "track_origin_y") origin.SetY(v);
    else if (n == "emission_z") emission_z = v;
    else if (n == "detector_theta") g.detector_theta = v;
    else if (n == "detector_phi") g.detector_phi = v;
    else if (n == "detector_delta_theta") g.detector_theta += v;
    else if (n == "detector_delta_phi") g.detector_phi += v;
    else if (n == "assembly_rotation_x") g.assembly_rotation_x = v;
    else if (n == "assembly_rotation_y") g.assembly_rotation_y = v;
    else if (n == "assembly_rotation_x_pivot_y") g.assembly_rotation_x_pivot.SetY(v);
    else if (n == "assembly_rotation_x_pivot_z") g.assembly_rotation_x_pivot.SetZ(v);
    else if (n == "mirror_delta_x") mirror_delta.SetX(v);
    else if (n == "mirror_delta_y") mirror_delta.SetY(v);
    else if (n == "mirror_delta_z") mirror_delta.SetZ(v);
    else if (n == "mirror_rotation_x") mirror_rotation.SetX(v);
    else if (n == "mirror_rotation_y") mirror_rotation.SetY(v);
    else if (n == "mirror_rotation_z") mirror_rotation.SetZ(v);
    else if (n == "detector_delta_x") detector_delta.SetX(v);
    else if (n == "detector_delta_y") detector_delta.SetY(v);
    else if (n == "detector_delta_z") detector_delta.SetZ(v);
    else if (n == "detector_rotation_x") detector_rotation.SetX(v);
    else if (n == "detector_rotation_y") detector_rotation.SetY(v);
    else if (n == "detector_rotation_z") detector_rotation.SetZ(v);
  }
  g.mirror_center += mirror_delta;
  g.detector_center += detector_delta;
  g.mirror_rotation_vector += mirror_rotation;
  g.detector_rotation_vector += detector_rotation;
  g.emission_z = emission_z;
  g.track_origin = origin;
  irt::update_track(g, eta, phi);
  irt::apply_inverse_assembly_to_track(g);
}

bool ring_hit(const ring_t &ring, const hit_t &hit)
{
  if (!std::isfinite(hit.x) || !std::isfinite(hit.y) ||
      !std::isfinite(hit.time) || ring.radius <= 0.)
    return false;
  const double dx = hit.x - ring.x0;
  const double dy = hit.y - ring.y0;
  const double c = std::cos(ring.phi);
  const double s = std::sin(ring.phi);
  const double minor = ring.radius *
                       std::sqrt(std::max(1.e-12,
                                          1. - ring.eccentricity * ring.eccentricity));
  const double xp = c * dx + s * dy;
  const double yp = -s * dx + c * dy;
  const double residual =
      std::abs(std::hypot(xp / ring.radius, yp / minor) - 1.) * ring.radius;
  return residual <= 4. * 1.5 &&
         std::abs(hit.time - ring.time) <= 4. * 1.;
}

} // namespace

int main(int argc, char **argv)
{
  std::string input;
  std::string output;
  std::string config;
  std::string fit_config;
  Long64_t max_events = -1;
  Long64_t fit_max_hits = 2000;
  bool fit_geometry = false;
  bool fit_track = false;
  bool fit_track_phi = false;
  bool fit_detector = false;
  bool fit_transverse = false;
  bool fit_harmonics = false;
  bool fit_mirror_rotation = false;
  bool fit_mirror_track = false;
  bool fit_mirror_phi = false;
  bool harmonic_objective = false;
  bool diagnostic = false;
  double target_theta = 0.037921467;
  int require_rings = -1;
  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help")
    ("input", po::value<std::string>(&input)->required(), "triggered ROOT input")
    ("output", po::value<std::string>(&output)->required(), "ROOT output")
    ("config", po::value<std::string>(&config)->required(), "nominal geometry configuration")
    ("fit-config", po::value<std::string>(&fit_config), "Minuit fit parameter configuration")
    ("max-events", po::value<Long64_t>(&max_events)->default_value(-1),
     "maximum frames to process")
    ("fit-max-hits", po::value<Long64_t>(&fit_max_hits)->default_value(fit_max_hits),
     "maximum deterministic Cherenkov hits used by the geometry fit")
    ("fit", po::bool_switch(&fit_geometry), "fit nominal geometry before filling output")
    ("fit-track", po::bool_switch(&fit_track), "fit only track eta and phi")
    ("fit-track-phi", po::bool_switch(&fit_track_phi), "fit only track phi, keeping configured eta fixed")
    ("fit-detector", po::bool_switch(&fit_detector), "fit detector theta and tilt x")
    ("fit-transverse", po::bool_switch(&fit_transverse), "fit mirror and detector transverse x/y shifts")
    ("fit-harmonics", po::bool_switch(&fit_harmonics), "fit detector orientation to minimize theta-phi harmonics")
    ("fit-mirror-rotation", po::bool_switch(&fit_mirror_rotation), "fit mirror rotation angle about its configured pivot")
    ("fit-mirror-track", po::bool_switch(&fit_mirror_track), "fit mirror rotation vector and track eta/phi together")
    ("fit-mirror-phi", po::bool_switch(&fit_mirror_phi), "fit mirror rotation vector and track phi with eta fixed to 3.16")
    ("harmonic-objective", po::bool_switch(&harmonic_objective), "include theta-versus-phi harmonic power in configured fit")
    ("diagnostic", po::bool_switch(&diagnostic), "scan one geometry parameter at a time")
    ("expected-cherenkov-angle", po::value<double>(&target_theta)->default_value(target_theta),
     "expected Cherenkov angle used by the geometry fit, in radians")
    ("target-theta", po::value<double>(&target_theta),
     "alias for --expected-cherenkov-angle")
    ("require-rings", po::value<int>(&require_rings)->default_value(require_rings),
     "require exactly this many accepted rings per frame")
    ;
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) { std::cout << options << '\n'; return 0; }
    po::notify(vm);
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n' << options << '\n';
    return 2;
  }

  trigger_reader_t reader;
  if (!reader.open(input))
    return 1;
  if (require_rings >= 0 && !reader.has_rings()) {
    std::cerr << "ERROR: --require-rings needs a ring tree in the input\n";
    return 1;
  }
  std::vector<std::pair<double, double>> points;
  std::vector<std::pair<double, double>> fit_points;
  Long64_t frames_read = 0;
  while (reader.next_spill() && (max_events < 0 || frames_read < max_events)) {
    while (reader.next_frame() && (max_events < 0 || frames_read < max_events)) {
      ++frames_read;
      if (require_rings >= 0 &&
          static_cast<int>(reader.rings().size()) != require_rings)
        continue;
      for (const auto &hit : reader.cherenkov_hits()) {
        if (std::isfinite(hit.x) && std::isfinite(hit.y)) {
          points.emplace_back(hit.x, hit.y);
          if (reader.has_rings()) {
            for (const auto &ring : reader.rings()) {
              if (ring_hit(ring, hit)) {
                fit_points.emplace_back(hit.x, hit.y);
                break;
              }
            }
          }
        }
      }
    }
  }

  const auto nominal = irt::nominal_geometry();
  auto geometry = nominal;
  if (!read_geometry_config(config, geometry)) {
    std::cerr << "ERROR: could not read valid geometry configuration: " << config << '\n';
    return 1;
  }
  irt::update_track(geometry, geometry.track_eta, geometry.track_phi);
  irt::apply_inverse_assembly_to_track(geometry);
  geometry.mirror_center += TVector3(2.11642, -12.6748, 0.);
  geometry.detector_center += TVector3(0.197618, 0.00373014, 0.);
  irt::geometry_t fitted_geometry = geometry;
  if (!fit_config.empty() && !fit_points.empty() && fit_max_hits > 0) {
    std::vector<fit_parameter_t> parameters;
    if (!read_fit_config(fit_config, parameters)) {
      std::cerr << "ERROR: invalid fit configuration: " << fit_config << '\n';
      return 1;
    }
    const std::size_t fit_count = std::min<std::size_t>(fit_points.size(), fit_max_hits);
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      apply_fit_parameters(candidate, parameters, par);
      double sum = 0.; int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first, fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.; ++used;
      }
      if (used == 0) return 1.e12;
      if (!harmonic_objective) return sum;
      TMatrixD normal(5, 5); TVectorD rhs(5);
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first, fit_points[i].second);
        if (!photon.valid) continue;
        const double b[5] = {1., std::cos(photon.phi), std::sin(photon.phi),
                             std::cos(2. * photon.phi), std::sin(2. * photon.phi)};
        for (int j = 0; j < 5; ++j) {
          rhs[j] += b[j] * (photon.theta - target_theta);
          for (int k = 0; k < 5; ++k) normal(j, k) += b[j] * b[k];
        }
      }
      const TVectorD h = normal.Invert() * rhs;
      return sum + 1.e8 * (h[1]*h[1] + h[2]*h[2] + h[3]*h[3] + h[4]*h[4]);
    };
    std::vector<double> start; start.reserve(parameters.size());
    for (const auto &p : parameters) start.push_back(p.start);
    ROOT::Math::Functor function(objective, static_cast<unsigned int>(parameters.size()));
    ROOT::Fit::Fitter fitter; fitter.SetFCN(function, start.data());
    for (std::size_t i = 0; i < parameters.size(); ++i) {
      fitter.Config().ParSettings(i).SetLimits(parameters[i].min, parameters[i].max);
      if (!parameters[i].free) fitter.Config().ParSettings(i).Fix();
    }
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: configured geometry fit did not converge\n";
      fitter.Result().Print(std::cerr);
    }
    else {
      const auto result = fitter.Result(); result.Print(std::cout);
      apply_fit_parameters(fitted_geometry, parameters, result.Parameters().data());
    }
  }
  auto calculate_harmonics = [&](const irt::geometry_t &candidate) {
    TMatrixD normal(5, 5);
    TVectorD rhs(5);
    for (const auto &point : fit_points) {
      const auto photon = irt::reconstruct(candidate, point.first, point.second);
      if (!photon.valid) continue;
      const double basis[5] = {1., std::cos(photon.phi), std::sin(photon.phi),
                               std::cos(2. * photon.phi), std::sin(2. * photon.phi)};
      for (int i = 0; i < 5; ++i) {
        rhs[i] += basis[i] * (photon.theta - target_theta);
        for (int j = 0; j < 5; ++j) normal(i, j) += basis[i] * basis[j];
      }
    }
    return normal.Invert() * rhs;
  };
  if (fit_track && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      irt::update_track(candidate, par[0], par[1]);
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[2] = {3.09673, .262884};
    ROOT::Math::Functor function(objective, 2);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    fitter.Config().ParSettings(0).SetLimits(2.5, 4.0);
    fitter.Config().ParSettings(1).SetLimits(-.5, .5);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: track fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      irt::update_track(fitted_geometry, result.Parameter(0), result.Parameter(1));
      std::cout << "fitted track_eta=" << result.Parameter(0)
                << " track_phi=" << result.Parameter(1) << '\n';
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_track_phi && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      irt::update_track(candidate, candidate.track_eta, par[0]);
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[1] = {geometry.track_phi};
    ROOT::Math::Functor function(objective, 1);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    fitter.Config().ParSettings(0).SetLimits(-.5, .5);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: track-phi fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      irt::update_track(fitted_geometry, geometry.track_eta, result.Parameter(0));
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fixed track_eta=" << geometry.track_eta
                << "; fitted track_phi=" << result.Parameter(0) << '\n'
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_detector && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.detector_theta = par[0];
      candidate.detector_tilt_x = par[1];
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[2] = {geometry.detector_theta, 0.};
    ROOT::Math::Functor function(objective, 2);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    fitter.Config().ParSettings(0).SetLimits(-.8, -.5);
    fitter.Config().ParSettings(1).SetLimits(-.02, .02);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: detector fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.detector_theta = result.Parameter(0);
      fitted_geometry.detector_tilt_x = result.Parameter(1);
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fitted detector_theta=" << result.Parameter(0)
                << " detector_tilt_x=" << result.Parameter(1) << '\n'
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_transverse && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.mirror_center += TVector3(par[0], par[1], 0.);
      candidate.detector_center += TVector3(par[2], par[3], 0.);
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[4] = {0., 0., 0., 0.};
    ROOT::Math::Functor function(objective, 4);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    for (int i = 0; i < 4; ++i)
      fitter.Config().ParSettings(i).SetLimits(-100., 100.);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: transverse fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.mirror_center += TVector3(result.Parameter(0), result.Parameter(1), 0.);
      fitted_geometry.detector_center += TVector3(result.Parameter(2), result.Parameter(3), 0.);
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fitted transverse shifts: mirror=("
                << result.Parameter(0) << "," << result.Parameter(1)
                << ") detector=(" << result.Parameter(2) << ","
                << result.Parameter(3) << ")\n"
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_harmonics && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto harmonic_objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.detector_theta = par[0];
      candidate.detector_phi = par[1];
      candidate.detector_tilt_x = par[2];
      candidate.detector_tilt_y = par[3];
      TMatrixD normal(5, 5);
      TVectorD rhs(5);
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double basis[5] = {1., std::cos(photon.phi), std::sin(photon.phi),
                                 std::cos(2. * photon.phi), std::sin(2. * photon.phi)};
        for (int j = 0; j < 5; ++j) {
          rhs[j] += basis[j] * (photon.theta - target_theta);
          for (int k = 0; k < 5; ++k) normal(j, k) += basis[j] * basis[k];
        }
      }
      const TVectorD h = normal.Invert() * rhs;
      return h[1] * h[1] + h[2] * h[2] + h[3] * h[3] + h[4] * h[4];
    };
    double start[4] = {geometry.detector_theta, geometry.detector_phi, 0., 0.};
    ROOT::Math::Functor function(harmonic_objective, 4);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    fitter.Config().ParSettings(0).SetLimits(-.8, -.5);
    fitter.Config().ParSettings(1).SetLimits(-.2, .2);
    fitter.Config().ParSettings(2).SetLimits(-.02, .02);
    fitter.Config().ParSettings(3).SetLimits(-.02, .02);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: harmonic fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.detector_theta = result.Parameter(0);
      fitted_geometry.detector_phi = result.Parameter(1);
      fitted_geometry.detector_tilt_x = result.Parameter(2);
      fitted_geometry.detector_tilt_y = result.Parameter(3);
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fitted harmonic objective=" << harmonic_objective(result.Parameters().data()) << '\n'
                << "fitted detector geometry: theta=" << result.Parameter(0)
                << " phi=" << result.Parameter(1)
                << " tilt_x=" << result.Parameter(2)
                << " tilt_y=" << result.Parameter(3) << '\n'
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_mirror_rotation && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.mirror_rotation_vector.SetXYZ(par[0], par[1], par[2]);
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[3] = {0., 0., 0.};
    ROOT::Math::Functor function(objective, 3);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    for (int i = 0; i < 3; ++i)
      fitter.Config().ParSettings(i).SetLimits(-.1, .1);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: mirror rotation fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.mirror_rotation_vector.SetXYZ(result.Parameter(0),
                                                    result.Parameter(1),
                                                    result.Parameter(2));
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fitted mirror rotation vector=(" << result.Parameter(0)
                << "," << result.Parameter(1) << "," << result.Parameter(2)
                << ") total angle=" << fitted_geometry.mirror_rotation_vector.Mag() << '\n'
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_mirror_track && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.mirror_rotation_vector.SetXYZ(par[0], par[1], par[2]);
      candidate.track.SetPtEtaPhi(1., par[3], par[4]);
      candidate.track = candidate.track.Unit();
      candidate.emission = 2639. * candidate.track;
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[5] = {0., 0., 0., 3.09673, .262884};
    ROOT::Math::Functor function(objective, 5);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    for (int i = 0; i < 3; ++i)
      fitter.Config().ParSettings(i).SetLimits(-.1, .1);
    fitter.Config().ParSettings(3).SetLimits(2.5, 4.0);
    fitter.Config().ParSettings(4).SetLimits(-.5, .5);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: mirror-track fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.mirror_rotation_vector.SetXYZ(result.Parameter(0),
                                                    result.Parameter(1),
                                                    result.Parameter(2));
      fitted_geometry.track.SetPtEtaPhi(1., result.Parameter(3), result.Parameter(4));
      fitted_geometry.track = fitted_geometry.track.Unit();
      fitted_geometry.emission = 2639. * fitted_geometry.track;
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fitted mirror rotation vector=("
                << result.Parameter(0) << "," << result.Parameter(1) << ","
                << result.Parameter(2) << ") track_eta=" << result.Parameter(3)
                << " track_phi=" << result.Parameter(4) << '\n'
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (fit_mirror_phi && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.mirror_rotation_vector.SetXYZ(par[0], par[1], par[2]);
      candidate.track.SetPtEtaPhi(1., 3.16, par[3]);
      candidate.track = candidate.track.Unit();
      candidate.emission = 2639. * candidate.track;
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto photon = irt::reconstruct(candidate, fit_points[i].first,
                                              fit_points[i].second);
        if (!photon.valid) continue;
        const double pull = (photon.theta - target_theta) / .001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[4] = {0., 0., 0., .262884};
    ROOT::Math::Functor function(objective, 4);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    for (int i = 0; i < 3; ++i)
      fitter.Config().ParSettings(i).SetLimits(-.1, .1);
    fitter.Config().ParSettings(3).SetLimits(-.5, .5);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: mirror-phi fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.mirror_rotation_vector.SetXYZ(result.Parameter(0),
                                                    result.Parameter(1),
                                                    result.Parameter(2));
      fitted_geometry.track.SetPtEtaPhi(1., 3.16, result.Parameter(3));
      fitted_geometry.track = fitted_geometry.track.Unit();
      fitted_geometry.emission = 2639. * fitted_geometry.track;
      const auto harmonics = calculate_harmonics(fitted_geometry);
      std::cout << "fixed track_eta=3.16; mirror rotation vector=("
                << result.Parameter(0) << "," << result.Parameter(1) << ","
                << result.Parameter(2) << ") track_phi=" << result.Parameter(3)
                << '\n'
                << "fitted theta-phi harmonics:"
                << " A0=" << harmonics[0]
                << " Ac1=" << harmonics[1]
                << " As1=" << harmonics[2]
                << " Ac2=" << harmonics[3]
                << " As2=" << harmonics[4] << '\n';
    }
  }
  if (diagnostic && !fit_points.empty()) {
    auto harmonics = [&](const irt::geometry_t &candidate) {
      TMatrixD normal(5, 5);
      TVectorD rhs(5);
      for (const auto &point : fit_points) {
        const auto photon = irt::reconstruct(candidate, point.first, point.second);
        if (!photon.valid) continue;
        const double phi = photon.phi;
        const double basis[5] = {1., std::cos(phi), std::sin(phi),
                                 std::cos(2. * phi), std::sin(2. * phi)};
        for (int i = 0; i < 5; ++i) {
          rhs[i] += basis[i] * (photon.theta - target_theta);
          for (int j = 0; j < 5; ++j)
            normal(i, j) += basis[i] * basis[j];
        }
      }
      TVectorD result(5);
      result = normal.Invert() * rhs;
      return result;
    };
    const auto nominal_harmonics = harmonics(geometry);
    std::cout << "nominal theta-phi harmonics:"
              << " A0=" << nominal_harmonics[0]
              << " Ac1=" << nominal_harmonics[1]
              << " As1=" << nominal_harmonics[2]
              << " Ac2=" << nominal_harmonics[3]
              << " As2=" << nominal_harmonics[4] << '\n';
    auto geometry_at = [&](int parameter, double value) {
      auto candidate = geometry;
      if (parameter == 0) candidate.emission = value * candidate.track;
      if (parameter == 1) { candidate.track.SetPtEtaPhi(1., value, .263759); candidate.track = candidate.track.Unit(); candidate.emission = 2639. * candidate.track; }
      if (parameter == 2) { candidate.track.SetPtEtaPhi(1., 3.1, value); candidate.track = candidate.track.Unit(); candidate.emission = 2639. * candidate.track; }
      if (parameter == 3) candidate.mirror_center += TVector3(value, 0., 0.);
      if (parameter == 4) candidate.mirror_center += TVector3(0., value, 0.);
      if (parameter == 5) candidate.mirror_center += TVector3(0., 0., value);
      if (parameter == 6) candidate.detector_center += TVector3(0., 0., value);
      if (parameter == 7) candidate.detector_center += TVector3(value, 0., 0.);
      if (parameter == 8) candidate.detector_center += TVector3(0., value, 0.);
      if (parameter == 9) candidate.detector_theta = value;
      if (parameter == 10) candidate.detector_phi = value;
      if (parameter == 11) candidate.detector_tilt_x = value;
      if (parameter == 12) candidate.detector_tilt_y = value;
      return candidate;
    };
    struct scan_t { const char *name; double lo; double hi; int parameter; };
    const scan_t scans[] = {
      {"emission_z", 2400., 2800., 0},
      {"track_eta", 2.5, 4.0, 1},
      {"track_phi", -.5, .5, 2},
      {"mirror_delta_x", -100., 100., 3},
      {"mirror_delta_y", -100., 100., 4},
      {"mirror_delta_z", -100., 100., 5},
      {"detector_center_delta_z", -100., 100., 6},
      {"detector_center_delta_x", -100., 100., 7},
      {"detector_center_delta_y", -100., 100., 8},
      {"detector_theta", -.8, -.5, 9},
      {"detector_phi", -.2, .2, 10},
      {"detector_tilt_x", -.02, .02, 11},
      {"detector_tilt_y", -.02, .02, 12}
    };
    for (const auto &scan : scans) {
      double best_value = scan.lo, best_objective = 1.e300;
      std::vector<double> scan_objectives(81, 1.e300);
      for (int istep = 0; istep <= 80; ++istep) {
        const double value = scan.lo + (scan.hi - scan.lo) * istep / 80.;
        auto candidate = geometry_at(scan.parameter, value);
        double objective = 0.;
        int used = 0;
        for (const auto &point : fit_points) {
          const auto photon = irt::reconstruct(candidate, point.first, point.second);
          if (!photon.valid) continue;
          const double pull = (photon.theta - target_theta) / .001;
          const double a = std::abs(pull);
          objective += a <= 3. ? pull * pull : 6. * a - 9.;
          ++used;
        }
        scan_objectives[istep] = used > 0 ? objective : 1.e300;
        if (used > 0 && objective < best_objective) {
          best_objective = objective;
          best_value = value;
        }
      }
      const double step = (scan.hi - scan.lo) / 80.;
      const double derivative =
          (scan_objectives[41] - scan_objectives[39]) / (2. * step);
      std::cout << "diagnostic " << scan.name
                << ": best=" << best_value
                << " objective=" << best_objective
                << " dobjective/dparameter=" << derivative
                << " |derivative|=" << std::abs(derivative) << '\n';
    }
  }
  if (fit_geometry && !fit_points.empty() && fit_max_hits > 0) {
    const std::size_t fit_count = std::min<std::size_t>(
        fit_points.size(), static_cast<std::size_t>(fit_max_hits));
    auto objective = [&](const double *par) {
      auto candidate = geometry;
      candidate.detector_tilt_x += par[0];
      candidate.detector_tilt_y += par[1];
      double sum = 0.;
      int used = 0;
      for (std::size_t i = 0; i < fit_count; ++i) {
        const auto &point = fit_points[i];
        const auto photon = irt::reconstruct(candidate, point.first, point.second);
        if (!photon.valid)
          continue;
        const double pull = (photon.theta - target_theta) / 0.001;
        const double a = std::abs(pull);
        sum += a <= 3. ? pull * pull : 6. * a - 9.;
        ++used;
      }
      return used > 0 ? sum : 1.e12;
    };
    double start[2] = {0., 0.};
    ROOT::Math::Functor function(objective, 2);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(function, start);
    fitter.Config().ParSettings(0).SetLimits(-100., 100.);
    fitter.Config().ParSettings(1).SetLimits(-100., 100.);
    if (!fitter.FitFCN()) {
      std::cerr << "WARNING: nominal geometry fit did not converge\n";
    } else {
      const auto result = fitter.Result();
      result.Print(std::cout);
      fitted_geometry.detector_tilt_x += result.Parameter(0);
      fitted_geometry.detector_tilt_y += result.Parameter(1);
    }
  }

  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << output << '\n';
    return 1;
  }

  auto h_angle = new TH1D("hCherenkovTheta", "reconstructed Cherenkov angle;#theta [rad];hits", 1000, 0., .1);
  auto h_phi = new TH1D("hCherenkovPhi", "reconstructed Cherenkov azimuth;phi [rad];hits", 720, -M_PI, M_PI);
  auto h_theta_phi = new TH2D("hCherenkovThetaPhi", "reconstructed Cherenkov angles;#phi [rad];#theta [rad]", 720, -M_PI, M_PI, 1000, 0., .1);
  auto h_reflection = new TH2D("hMirrorXY", "mirror reflection points;x [mm];y [mm]", 1000, -1000., 1000., 1000, -1000., 1000.);
  auto h_detector = new TH2D("hDetectorXY", "detector hits;x [mm];y [mm]", 400, -200., 200., 400, -200., 200.);

  // Recreate all frame-aligned trees entry-by-entry, appending the IRT angles
  // to the Cherenkov tree while preserving the original branch contents.
  std::unique_ptr<TFile> input_file(TFile::Open(input.c_str(), "READ"));
  TTree *input_cherenkov = input_file ?
      dynamic_cast<TTree *>(input_file->Get("cherenkov")) : nullptr;
  TTree *output_cherenkov = nullptr;
  UShort_t cherenkov_nhits = 0;
  Float_t cherenkov_x[65535] = {};
  Float_t cherenkov_y[65535] = {};
  Float_t cherenkov_theta[65535] = {};
  Float_t cherenkov_phi[65535] = {};
  if (!input_cherenkov ||
      input_cherenkov->SetBranchAddress("nhits", &cherenkov_nhits) < 0 ||
      input_cherenkov->SetBranchAddress("x", cherenkov_x) < 0 ||
      input_cherenkov->SetBranchAddress("y", cherenkov_y) < 0) {
    std::cerr << "ERROR: input does not contain a usable cherenkov tree\n";
    return 1;
  }
  file->cd();
  output_cherenkov = input_cherenkov->CloneTree(0);
  output_cherenkov->SetName("cherenkov");
  auto output_irt = new TTree("irt", "IRT reconstruction, one entry per frame");
  output_irt->Branch("nhits", &cherenkov_nhits, "nhits/s");
  output_irt->Branch("theta", cherenkov_theta, "theta[nhits]/F");
  output_irt->Branch("phi", cherenkov_phi, "phi[nhits]/F");

  TTree *input_frames = dynamic_cast<TTree *>(input_file->Get("frames"));
  TTree *input_ring = dynamic_cast<TTree *>(input_file->Get("ring"));
  TTree *output_frames = input_frames ? input_frames->CloneTree(0) : nullptr;
  TTree *output_ring = input_ring ? input_ring->CloneTree(0) : nullptr;

  TIter key_iterator(input_file->GetListOfKeys());
  while (auto *key = dynamic_cast<TKey *>(key_iterator())) {
    const std::string name = key->GetName();
    if (name == "frames" || name == "ring" || name == "cherenkov")
      continue;
    auto *tree = dynamic_cast<TTree *>(key->ReadObj());
    if (!tree)
      continue;
    file->cd();
    auto *copy = tree->CloneTree(-1, "fast");
    if (!copy) {
      std::cerr << "ERROR: could not clone input tree '" << name << "'\n";
      return 1;
    }
    copy->Write();
  }

  Long64_t frames = 0, hits = 0, valid = 0;
  for (const auto &point : points) {
        ++hits;
        const auto photon = irt::reconstruct(fitted_geometry, point.first, point.second);
        if (!photon.valid)
          continue;
        ++valid;
        h_angle->Fill(photon.theta);
        h_phi->Fill(photon.phi);
        h_theta_phi->Fill(photon.phi, photon.theta);
        h_reflection->Fill(photon.reflection.X(), photon.reflection.Y());
        h_detector->Fill(point.first, point.second);
  }
  frames = frames_read;
  for (Long64_t entry = 0; entry < input_cherenkov->GetEntries(); ++entry) {
    if (input_frames && input_frames->GetEntry(entry) <= 0)
      return 1;
    if (input_ring && input_ring->GetEntry(entry) <= 0)
      return 1;
    if (input_cherenkov->GetEntry(entry) <= 0)
      return 1;
    if (cherenkov_nhits > 65535) {
      std::cerr << "ERROR: cherenkov hit count exceeds output capacity\n";
      return 1;
    }
    for (unsigned int i = 0; i < cherenkov_nhits; ++i) {
      const auto photon = irt::reconstruct(fitted_geometry,
                                            cherenkov_x[i], cherenkov_y[i]);
      cherenkov_theta[i] = photon.valid ? static_cast<Float_t>(photon.theta) :
                                           std::numeric_limits<Float_t>::quiet_NaN();
      cherenkov_phi[i] = photon.valid ? static_cast<Float_t>(photon.phi) :
                                         std::numeric_limits<Float_t>::quiet_NaN();
    }
    if (output_frames && output_frames->Fill() < 0)
      return 1;
    if (output_ring && output_ring->Fill() < 0)
      return 1;
    output_cherenkov->Fill();
    output_irt->Fill();
  }
  file->Write("", TObject::kOverwrite);
  h_angle->Write();
  h_phi->Write();
  h_theta_phi->Write();
  h_reflection->Write();
  h_detector->Write();
  file->Close();
  std::cout << "frames processed: " << frames << '\n'
            << "Cherenkov hits:  " << hits << '\n'
            << "valid photons:   " << valid << '\n'
            << "output:           " << output << '\n';
  return 0;
}

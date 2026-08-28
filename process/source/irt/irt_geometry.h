#pragma once

#include <TVector3.h>

#include <cmath>

namespace irt {

struct geometry_t {
  TVector3 track;
  TVector3 track_origin;
  TVector3 emission;
  double track_eta = 0.;
  double track_phi = 0.;
  double emission_z = 0.;
  TVector3 mirror_center;
  double mirror_radius = 0.;
  TVector3 mirror_pivot;
  TVector3 mirror_rotation_vector;
  TVector3 detector_center;
  double detector_radius = 0.;
  double detector_theta = 0.;
  double detector_phi = 0.;
  double detector_tilt_x = 0.;
  double detector_tilt_y = 0.;
  TVector3 detector_rotation_vector;
  TVector3 assembly_rotation_x_pivot;
  TVector3 assembly_rotation_y_pivot;
  double assembly_rotation_x = 0.;
  double assembly_rotation_y = 0.;
};

struct photon_t {
  bool valid = false;
  TVector3 detector;
  TVector3 reflection;
  double theta = 0.;
  double phi = 0.;
};

inline bool finite(const TVector3 &v)
{
  return std::isfinite(v.X()) && std::isfinite(v.Y()) &&
         std::isfinite(v.Z());
}

inline bool inverse_ray_trace(const TVector3 &emission,
                              const TVector3 &detector,
                              const TVector3 &center,
                              double radius,
                              TVector3 &reflection)
{
  const TVector3 ce = emission - center;
  const TVector3 cd = detector - center;
  const double a = ce.Mag();
  const double d = cd.Mag();
  const double alpha = ce.Angle(cd);
  const double sine = std::sin(alpha);
  if (!(a > 0.) || !(d > 0.) || !(radius > 0.) ||
      !std::isfinite(alpha) || std::abs(sine) < 1.e-10)
    return false;

  auto function = [&](double beta) {
    return a * d * std::sin(alpha - 2. * beta) +
           radius * (a * std::sin(beta) -
                     d * std::sin(alpha - beta));
  };
  auto derivative = [&](double beta) {
    return -2. * a * d * std::cos(alpha - 2. * beta) +
           radius * (a * std::cos(beta) + d * std::cos(alpha - beta));
  };

  double beta = 0.5 * alpha;
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double derivative_value = derivative(beta);
    if (!std::isfinite(derivative_value) ||
        std::abs(derivative_value) < 1.e-12)
      return false;
    const double delta = function(beta) / derivative_value;
    if (!std::isfinite(delta) || std::abs(delta) > 1.)
      return false;
    beta -= delta;
    if (std::abs(delta) < 1.e-10)
      break;
  }

  reflection = center +
    (radius * std::cos(beta) / a -
     radius * std::sin(beta) * std::cos(alpha) / (a * sine)) * ce +
    radius * std::sin(beta) / (d * sine) * cd;
  return finite(reflection);
}

inline photon_t reconstruct(const geometry_t &geometry, double x, double y)
{
  photon_t result;
  TVector3 normal;
  normal.SetMagThetaPhi(1., geometry.detector_theta, geometry.detector_phi);
  TVector3 plane_center = geometry.detector_center + geometry.detector_radius * normal;
  if (geometry.detector_rotation_vector.Mag2() > 0.)
    normal.Rotate(geometry.detector_rotation_vector.Mag(),
                  geometry.detector_rotation_vector.Unit());
  TVector3 detector_reference(0., 1., 0.);
  if (std::abs(detector_reference.Dot(normal)) > .999)
    detector_reference.SetXYZ(1., 0., 0.);
  const TVector3 detector_axis_x =
      (detector_reference - detector_reference.Dot(normal) * normal).Unit();
  const TVector3 detector_axis_y = normal.Cross(detector_axis_x).Unit();
  normal.Rotate(geometry.detector_tilt_x, detector_axis_x);
  normal.Rotate(geometry.detector_tilt_y, detector_axis_y);
  result.detector = plane_center;
  TVector3 local(x, y, 0.);
  local.RotateUz(normal);
  result.detector += local;
  TVector3 mirror_center = geometry.mirror_center;
  if (geometry.mirror_rotation_vector.Mag2() > 0.) {
    mirror_center -= geometry.mirror_pivot;
    mirror_center.Rotate(geometry.mirror_rotation_vector.Mag(),
                         geometry.mirror_rotation_vector.Unit());
    mirror_center += geometry.mirror_pivot;
  }
  if (!inverse_ray_trace(geometry.emission, result.detector,
                         mirror_center, geometry.mirror_radius,
                         result.reflection))
    return result;

  const TVector3 direction = (result.reflection - geometry.emission).Unit();
  if (!finite(direction) || direction.Dot(geometry.track) <= 0.)
    return result;
  result.theta = geometry.track.Angle(direction);
  TVector3 track_reference(0., 1., 0.);
  if (std::abs(track_reference.Dot(geometry.track)) > .999)
    track_reference.SetXYZ(1., 0., 0.);
  const TVector3 track_axis_x =
    (track_reference - track_reference.Dot(geometry.track) * geometry.track).Unit();
  const TVector3 track_axis_y = geometry.track.Cross(track_axis_x).Unit();
  result.phi = std::atan2(direction.Dot(track_axis_y), direction.Dot(track_axis_x));
  result.valid = std::isfinite(result.theta) && std::isfinite(result.phi);
  return result;
}

inline void update_track(geometry_t &geometry, double eta, double phi)
{
  geometry.track_eta = eta;
  geometry.track_phi = phi;
  geometry.track.SetPtEtaPhi(1., eta, phi);
  geometry.track = geometry.track.Unit();
  if (std::abs(geometry.track.Z()) > 1.e-12)
    geometry.emission = geometry.track_origin +
      ((geometry.emission_z - geometry.track_origin.Z()) / geometry.track.Z()) *
      geometry.track;
}

inline void apply_inverse_assembly_to_track(geometry_t &geometry)
{
  TVector3 direction(0., 0., 1.);
  TVector3 origin = geometry.track_origin;
  origin -= geometry.assembly_rotation_y_pivot;
  origin.Rotate(-geometry.assembly_rotation_y, TVector3(0., 1., 0.));
  origin += geometry.assembly_rotation_y_pivot;
  origin -= geometry.assembly_rotation_x_pivot;
  origin.Rotate(-geometry.assembly_rotation_x, TVector3(1., 0., 0.));
  origin += geometry.assembly_rotation_x_pivot;
  direction.Rotate(-geometry.assembly_rotation_y, TVector3(0., 1., 0.));
  direction.Rotate(-geometry.assembly_rotation_x, TVector3(1., 0., 0.));
  geometry.track_origin = origin;
  geometry.track = direction.Unit();
  geometry.emission = geometry.track_origin +
    ((geometry.emission_z - geometry.track_origin.Z()) / geometry.track.Z()) *
    geometry.track;
}

inline geometry_t nominal_geometry()
{
  geometry_t result;
  result.track.SetPtEtaPhi(1., 3.1, .263759);
  result.track = result.track.Unit();
  result.track_origin.SetXYZ(0., 0., 0.);
  result.emission_z = 2639.;
  update_track(result, 3.1, .263759);
  result.mirror_center.SetXYZ(1150.08, 0., 939.);
  result.mirror_radius = 2203.01;
  result.detector_center.SetXYZ(1838., 0., 1414.);
  result.detector_radius = 1100.;
  result.detector_theta = -.64805421;
  result.detector_phi = 0.;
  return result;
}

} // namespace irt

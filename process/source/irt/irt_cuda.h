#pragma once

namespace irt::cuda {

struct geometry_t {
  double track[3];
  double emission[3];
  double mirror_center[3];
  double mirror_pivot[3];
  double mirror_rotation[3];
  double mirror_radius;
  double detector_center[3];
  double detector_radius;
  double detector_normal[3];
  double detector_rotation[3];
  double detector_tilt_x;
  double detector_tilt_y;
};

bool available();

bool reconstruct(const geometry_t &geometry,
                 const float *x,
                 const float *y,
                 const float *time,
                 int count,
                 float *theta,
                 float *phi,
                 float *emission_time,
                 unsigned char *valid);

} // namespace irt::cuda

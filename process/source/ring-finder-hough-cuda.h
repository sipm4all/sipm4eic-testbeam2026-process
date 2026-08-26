#pragma once

#include <string>

namespace ring_hough_cuda {

struct grid_t {
  double min_x0;
  double min_y0;
  double min_radius;
  double min_t;
  double min_e;
  double min_phi;
  double x0_step;
  double y0_step;
  double radius_step;
  double t_step;
  double e_step;
  double phi_step;
  double spatial_resolution;
  double time_resolution;
  int nx;
  int ny;
  int nr;
  int nt;
  int ne;
  int nphi;
};

struct candidate_t {
  int x0_bin = -1;
  int y0_bin = -1;
  int radius_bin = -1;
  int t_bin = -1;
  int e_bin = -1;
  int phi_bin = -1;
  float score = 0.f;
};

class engine_t {
public:
  engine_t();
  ~engine_t();

  engine_t(const engine_t &) = delete;
  engine_t &operator=(const engine_t &) = delete;

  bool initialize(const grid_t &grid, std::string &error);
  bool find(const float *x, const float *y, const float *time,
            int nhits, int max_candidates,
            candidate_t *candidates, int &candidate_count,
            std::string &error);

private:
  void *impl_ = nullptr;
};

#if !defined(RING_HOUGH_HAS_CUDA) && \
    !defined(RING_HOUGH_CUDA_IMPLEMENTATION)
inline engine_t::engine_t() = default;
inline engine_t::~engine_t() = default;

inline bool
engine_t::initialize(const grid_t &, std::string &error)
{
  error = "CUDA backend is not compiled";
  return false;
}

inline bool
engine_t::find(const float *, const float *, const float *, int, int,
               candidate_t *, int &, std::string &error)
{
  error = "CUDA backend is not compiled";
  return false;
}
#endif

}

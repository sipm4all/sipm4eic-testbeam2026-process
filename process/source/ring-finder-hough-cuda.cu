#define RING_HOUGH_CUDA_IMPLEMENTATION
#include "ring-finder-hough-cuda.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <vector>

namespace {

constexpr int threads_per_block = 256;
// The full default time range can require about 1.1 GB at a 1 mm spatial
// grid. Keep an explicit upper bound while allowing that valid configuration.
constexpr std::size_t max_accumulator_entries = 300000000;
constexpr float gaussian_cut = 4.f;

struct device_state_t {
  ring_hough_cuda::grid_t grid{};
  float *x0_map = nullptr;
  float *y0_map = nullptr;
  float *radius_map = nullptr;
  float *time_map = nullptr;
  float *x = nullptr;
  float *y = nullptr;
  float *time = nullptr;
  float *accumulator = nullptr;
  unsigned long long *reduction_a = nullptr;
  int hit_capacity = 0;
  std::size_t accumulator_entries = 0;
  std::size_t reduction_capacity = 0;
};

struct ransac_device_state_t {
  float *x = nullptr;
  float *y = nullptr;
  float *time = nullptr;
  ring_hough_cuda::ransac_model_t *models = nullptr;
  ring_hough_cuda::ransac_stats_t *stats = nullptr;
  std::uint32_t *masks = nullptr;
  int hit_capacity = 0;
  int model_capacity = 0;
  std::size_t mask_capacity = 0;
};

std::string
cuda_error(const char *operation, cudaError_t error)
{
  std::ostringstream message;
  message << operation << ": " << cudaGetErrorString(error);
  return message.str();
}

__global__ void
coordinate_map_kernel(ring_hough_cuda::grid_t grid,
                      float *x0_map, float *y0_map,
                      float *radius_map, float *time_map)
{
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) *
                            blockDim.x + threadIdx.x;
  const std::size_t entries = static_cast<std::size_t>(grid.nx) * grid.ny *
                              grid.nr * grid.nt * grid.ne * grid.nphi;
  if (index >= entries)
    return;

  const int plane = grid.nr * grid.nt * grid.ne * grid.nphi;
  const int center = static_cast<int>(index / plane);
  const int remainder = static_cast<int>(index % plane);
  const int ix = center / grid.ny;
  const int iy = center % grid.ny;
  const int radius_plane = grid.nt * grid.ne * grid.nphi;
  const int ir = remainder / radius_plane;
  const int radius_remainder = remainder % radius_plane;
  const int it = radius_remainder / (grid.ne * grid.nphi);

  x0_map[index] = static_cast<float>(grid.min_x0 + ix * grid.x0_step);
  y0_map[index] = static_cast<float>(grid.min_y0 + iy * grid.y0_step);
  radius_map[index] = static_cast<float>(grid.min_radius +
                                         ir * grid.radius_step);
  time_map[index] = static_cast<float>(grid.min_t + it * grid.t_step);
}

__device__ float
circle_radial_residual(float x0, float y0, float radius, float x, float y)
{
  const float dx = x - x0;
  const float dy = y - y0;
  return hypotf(dx, dy) - radius;
}

__device__ float
ellipse_radial_residual(float x0, float y0, float radius, float eccentricity,
                        float phi, float x, float y)
{
  const float dx = x - x0;
  const float dy = y - y0;
  const float cosine = cosf(phi);
  const float sine = sinf(phi);
  const float xr = cosine * dx + sine * dy;
  const float yr = -sine * dx + cosine * dy;
  const float minor = radius * sqrtf(fmaxf(0.f, 1.f - eccentricity * eccentricity));
  const float angle = atan2f(yr, xr);
  const float ca = cosf(angle);
  const float sa = sinf(angle);
  const float denominator = hypotf(minor * ca, radius * sa);
  if (!(denominator > 0.f))
    return __int_as_float(0x7f800000);
  const float ellipse_radius = radius * minor / denominator;
  return hypotf(dx, dy) - ellipse_radius;
}

__global__ void
score_kernel(const float *x0_map, const float *y0_map,
             const float *radius_map, const float *time_map,
             const float *x, const float *y, const float *time, int nhits,
             ring_hough_cuda::grid_t grid, float *accumulator)
{
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) *
                            blockDim.x + threadIdx.x;
  const std::size_t entries = static_cast<std::size_t>(grid.nx) * grid.ny *
                              grid.nr * grid.nt * grid.ne * grid.nphi;
  if (index >= entries)
    return;

  const float x0 = x0_map[index];
  const float y0 = y0_map[index];
  const float radius = radius_map[index];
  const float ring_time = time_map[index];
  const float spatial_sigma = static_cast<float>(grid.spatial_resolution);
  const float time_sigma = static_cast<float>(grid.time_resolution);
  const int plane = grid.nr * grid.nt * grid.ne * grid.nphi;
  const int remainder = static_cast<int>(index % plane);
  const int radius_plane = grid.nt * grid.ne * grid.nphi;
  const int radius_remainder = remainder % radius_plane;
  const int ellipse_remainder = radius_remainder % (grid.ne * grid.nphi);
  const int eccentricity_bin = ellipse_remainder / grid.nphi;
  const int phi_bin = ellipse_remainder % grid.nphi;
  const float eccentricity = static_cast<float>(
      grid.min_e + eccentricity_bin * grid.e_step);
  const float phi = static_cast<float>(grid.min_phi + phi_bin * grid.phi_step);

  float score = 0.f;
  for (int hit = 0; hit < nhits; ++hit) {
    const float dr = ellipse_radial_residual(
        x0, y0, radius, eccentricity, phi, x[hit], y[hit]);
    const float dt = time[hit] - ring_time;
    if (fabsf(dr) > gaussian_cut * spatial_sigma ||
        fabsf(dt) > gaussian_cut * time_sigma)
      continue;
    const float spatial_pull = dr / spatial_sigma;
    const float time_pull = dt / time_sigma;
    score += expf(-0.5f * spatial_pull * spatial_pull) /
             (2.50662827463f * spatial_sigma) *
             expf(-0.5f * time_pull * time_pull) /
             (2.50662827463f * time_sigma);
  }
  accumulator[index] = score;
}

__global__ void
score_circle_kernel(const float *x0_map, const float *y0_map,
                    const float *radius_map, const float *time_map,
                    const float *x, const float *y, const float *time,
                    int nhits, ring_hough_cuda::grid_t grid,
                    float *accumulator)
{
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) *
                            blockDim.x + threadIdx.x;
  const std::size_t entries = static_cast<std::size_t>(grid.nx) * grid.ny *
                              grid.nr * grid.nt;
  if (index >= entries)
    return;

  const float x0 = x0_map[index];
  const float y0 = y0_map[index];
  const float radius = radius_map[index];
  const float ring_time = time_map[index];
  const float spatial_sigma = static_cast<float>(grid.spatial_resolution);
  const float time_sigma = static_cast<float>(grid.time_resolution);

  float score = 0.f;
  for (int hit = 0; hit < nhits; ++hit) {
    const float dr = circle_radial_residual(
        x0, y0, radius, x[hit], y[hit]);
    const float dt = time[hit] - ring_time;
    if (fabsf(dr) > gaussian_cut * spatial_sigma ||
        fabsf(dt) > gaussian_cut * time_sigma)
      continue;
    const float spatial_pull = dr / spatial_sigma;
    const float time_pull = dt / time_sigma;
    score += expf(-0.5f * spatial_pull * spatial_pull) /
             (2.50662827463f * spatial_sigma) *
             expf(-0.5f * time_pull * time_pull) /
             (2.50662827463f * time_sigma);
  }
  accumulator[index] = score;
}

__device__ bool
fit_circle_from_stats(const ring_hough_cuda::ransac_stats_t &stats,
                      ring_hough_cuda::ransac_model_t &model)
{
  double normal[3][3];
  double rhs[3];
  for (int row = 0; row < 3; ++row) {
    rhs[row] = stats.rhs[row];
    for (int column = 0; column < 3; ++column)
      normal[row][column] = stats.normal[row][column];
  }

  for (int column = 0; column < 3; ++column) {
    int pivot = column;
    for (int row = column + 1; row < 3; ++row) {
      if (fabs(normal[row][column]) > fabs(normal[pivot][column]))
        pivot = row;
    }
    if (fabs(normal[pivot][column]) < 1.e-12)
      return false;
    if (pivot != column) {
      for (int j = column; j < 3; ++j) {
        const double value = normal[column][j];
        normal[column][j] = normal[pivot][j];
        normal[pivot][j] = value;
      }
      const double value = rhs[column];
      rhs[column] = rhs[pivot];
      rhs[pivot] = value;
    }
    for (int row = column + 1; row < 3; ++row) {
      const double factor = normal[row][column] / normal[column][column];
      for (int j = column; j < 3; ++j)
        normal[row][j] -= factor * normal[column][j];
      rhs[row] -= factor * rhs[column];
    }
  }

  double solution[3] = {};
  for (int row = 2; row >= 0; --row) {
    double value = rhs[row];
    for (int column = row + 1; column < 3; ++column)
      value -= normal[row][column] * solution[column];
    solution[row] = value / normal[row][row];
  }
  model.x = static_cast<float>(0.5 * solution[0]);
  model.y = static_cast<float>(0.5 * solution[1]);
  const double radius_squared = solution[2] +
                                model.x * model.x + model.y * model.y;
  if (!(radius_squared > 0.) || !isfinite(radius_squared))
    return false;
  model.radius = static_cast<float>(sqrt(radius_squared));
  return isfinite(model.x) && isfinite(model.y) && isfinite(model.radius);
}

__global__ void
refine_ransac_kernel(
    const float *x, const float *y, const float *time, int nhits,
    ring_hough_cuda::ransac_model_t *models, int nmodels,
    double spatial_cut, double time_cut, bool write_masks, int mask_words,
    ring_hough_cuda::ransac_stats_t *stats, std::uint32_t *masks)
{
  constexpr int nterms = 14;
  const int model_index = static_cast<int>(blockIdx.x);
  if (model_index >= nmodels)
    return;

  __shared__ int partial_count[threads_per_block];
  __shared__ double partial[nterms][threads_per_block];
  __shared__ ring_hough_cuda::ransac_model_t model;
  __shared__ ring_hough_cuda::ransac_stats_t previous;
  const int thread = threadIdx.x;
  if (thread == 0) {
    model = models[model_index];
    previous = {};
  }
  __syncthreads();

  for (int phase = 0; phase < 5; ++phase) {
    if (thread == 0 && phase >= 2 && model.valid &&
        previous.inliers >= 3)
      model.ring_time = static_cast<float>(
          previous.sum_time / previous.inliers);
    __syncthreads();

    partial_count[thread] = 0;
    for (int term = 0; term < nterms; ++term)
      partial[term][thread] = 0.;
    __syncthreads();

    if (model.valid) {
      for (int hit = thread; hit < nhits; hit += blockDim.x) {
        const float dr = hypotf(x[hit] - model.x, y[hit] - model.y) -
                         model.radius;
        const float dt = time[hit] - model.ring_time;
        if (fabsf(dr) > spatial_cut || fabsf(dt) > time_cut)
          continue;

        const double hit_x = x[hit];
        const double hit_y = y[hit];
        const double value = hit_x * hit_x + hit_y * hit_y;
        ++partial_count[thread];
        partial[0][thread] += time[hit];
        partial[1][thread] += fabs(static_cast<double>(dr));
        partial[2][thread] += hit_x * hit_x;
        partial[3][thread] += hit_x * hit_y;
        partial[4][thread] += hit_x;
        partial[5][thread] += hit_x * hit_y;
        partial[6][thread] += hit_y * hit_y;
        partial[7][thread] += hit_y;
        partial[8][thread] += hit_x;
        partial[9][thread] += hit_y;
        partial[10][thread] += 1.;
        partial[11][thread] += hit_x * value;
        partial[12][thread] += hit_y * value;
        partial[13][thread] += value;
        if (write_masks && phase == 4)
          atomicOr(&masks[static_cast<std::size_t>(model_index) *
                          mask_words + hit / 32],
                   1u << (hit % 32));
      }
    }
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
      if (thread < stride) {
        partial_count[thread] += partial_count[thread + stride];
        for (int term = 0; term < nterms; ++term)
          partial[term][thread] += partial[term][thread + stride];
      }
      __syncthreads();
    }
    if (thread == 0) {
      ring_hough_cuda::ransac_stats_t current{};
      current.inliers = partial_count[0];
      current.sum_time = partial[0][0];
      current.residual_sum = partial[1][0];
      for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column)
          current.normal[row][column] =
              partial[2 + row * 3 + column][0];
        current.rhs[row] = partial[11 + row][0];
      }
      previous = current;
      if (current.inliers < 3)
        model.valid = 0;
      else if (phase == 0) {
        model.ring_time = static_cast<float>(
            current.sum_time / current.inliers);
      } else if (phase < 4 && !fit_circle_from_stats(current, model)) {
        model.valid = 0;
      }
      if (phase == 4)
        stats[model_index] = current;
    }
    __syncthreads();
  }
  if (thread == 0)
    models[model_index] = model;
}

__device__ unsigned long long
score_key(float score, std::size_t index)
{
  if (!(score > 0.f))
    return 0;
  // All scores are positive, so their IEEE-754 bit patterns have the same
  // ordering as the numerical values. The low word makes ties deterministic.
  const unsigned long long score_bits =
      static_cast<unsigned long long>(__float_as_uint(score));
  const unsigned long long inverse_index =
      0xffffffffull - static_cast<unsigned long long>(index);
  return (score_bits << 32) | inverse_index;
}

__global__ void
reduce_block_max_kernel(const float *scores, std::size_t entries,
                        unsigned long long *block_maxima)
{
  __shared__ unsigned long long keys[threads_per_block];
  const int thread = threadIdx.x;
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) *
                            blockDim.x + thread;
  keys[thread] = index < entries ? score_key(scores[index], index) : 0;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (thread < stride)
      keys[thread] = max(keys[thread], keys[thread + stride]);
    __syncthreads();
  }
  if (thread == 0)
    block_maxima[blockIdx.x] = keys[0];
}

__global__ void
suppress_neighborhood_kernel(float *accumulator,
                             ring_hough_cuda::grid_t grid,
                             int best_index)
{
  if (blockIdx.x != 0 || threadIdx.x != 0)
    return;

  const int plane = grid.nr * grid.nt * grid.ne * grid.nphi;
  const int center = best_index / plane;
  const int remainder = best_index % plane;
  const int ellipse_plane = grid.ne * grid.nphi;
  const int radius_plane = grid.nt * ellipse_plane;
  const int radius_bin = remainder / radius_plane;
  const int radius_remainder = remainder % radius_plane;
  const int time_bin = radius_remainder / ellipse_plane;
  const int ellipse_remainder = radius_remainder % ellipse_plane;
  const int eccentricity_bin = ellipse_remainder / grid.nphi;
  const int phi_bin = ellipse_remainder % grid.nphi;

  const int center_x = center / grid.ny;
  const int center_y = center % grid.ny;
  for (int dx = -1; dx <= 1; ++dx) {
    const int ix = center_x + dx;
    if (ix < 0 || ix >= grid.nx)
      continue;
    for (int dy = -1; dy <= 1; ++dy) {
      const int iy = center_y + dy;
      if (iy < 0 || iy >= grid.ny)
        continue;
      const int suppressed_center = ix * grid.ny + iy;
      for (int dr = -1; dr <= 1; ++dr) {
        const int ir = radius_bin + dr;
        if (ir < 0 || ir >= grid.nr)
          continue;
        for (int dt = -1; dt <= 1; ++dt) {
          const int it = time_bin + dt;
          if (it < 0 || it >= grid.nt)
            continue;
          for (int de = -1; de <= 1; ++de) {
            const int ie = eccentricity_bin + de;
            if (ie < 0 || ie >= grid.ne)
              continue;
            for (int dphi = -1; dphi <= 1; ++dphi) {
              const int iphi = phi_bin + dphi;
              if (iphi < 0 || iphi >= grid.nphi)
                continue;
              const std::size_t index =
                  ((static_cast<std::size_t>(suppressed_center) * grid.nr +
                    ir) * grid.nt + it) * grid.ne * grid.nphi +
                  static_cast<std::size_t>(ie) * grid.nphi + iphi;
              accumulator[index] = -1.f;
            }
          }
        }
      }
    }
  }
}

void
release(device_state_t &state)
{
  cudaFree(state.x0_map);
  cudaFree(state.y0_map);
  cudaFree(state.radius_map);
  cudaFree(state.time_map);
  cudaFree(state.x);
  cudaFree(state.y);
  cudaFree(state.time);
  cudaFree(state.accumulator);
  state.grid = {};
  state.x0_map = nullptr;
  state.y0_map = nullptr;
  state.radius_map = nullptr;
  state.time_map = nullptr;
  state.x = nullptr;
  state.y = nullptr;
  state.time = nullptr;
  state.accumulator = nullptr;
  cudaFree(state.reduction_a);
  state.reduction_a = nullptr;
  state.hit_capacity = 0;
  state.accumulator_entries = 0;
  state.reduction_capacity = 0;
}

void
release(ransac_device_state_t &state)
{
  cudaFree(state.x);
  cudaFree(state.y);
  cudaFree(state.time);
  cudaFree(state.models);
  cudaFree(state.stats);
  cudaFree(state.masks);
  state.x = nullptr;
  state.y = nullptr;
  state.time = nullptr;
  state.models = nullptr;
  state.stats = nullptr;
  state.masks = nullptr;
  state.hit_capacity = 0;
  state.model_capacity = 0;
  state.mask_capacity = 0;
}

bool
same_grid(const ring_hough_cuda::grid_t &first,
          const ring_hough_cuda::grid_t &second)
{
  return first.x0_step == second.x0_step &&
         first.y0_step == second.y0_step &&
         first.radius_step == second.radius_step &&
         first.t_step == second.t_step &&
         first.e_step == second.e_step &&
         first.phi_step == second.phi_step &&
         first.spatial_resolution == second.spatial_resolution &&
         first.time_resolution == second.time_resolution &&
         first.nx == second.nx && first.ny == second.ny &&
         first.nr == second.nr && first.nt == second.nt &&
         first.ne == second.ne && first.nphi == second.nphi &&
         first.circle_mode == second.circle_mode;
}

}

namespace ring_hough_cuda {

engine_t::engine_t() = default;

engine_t::~engine_t()
{
  auto *state = static_cast<device_state_t *>(impl_);
  if (state) {
    release(*state);
    delete state;
  }
}

bool
engine_t::initialize(const grid_t &grid, std::string &error)
{
  auto *old_state = static_cast<device_state_t *>(impl_);
  if (old_state) {
    if (same_grid(old_state->grid, grid)) {
      // Independent RANSAC seeds can have different local origins while
      // sharing the same accumulator shape. Keep the allocations and update
      // only the coordinate origin used by the kernels and result decoding.
      old_state->grid = grid;
      const std::size_t entries = old_state->accumulator_entries;
      const int blocks = static_cast<int>(
          (entries + threads_per_block - 1) / threads_per_block);
      coordinate_map_kernel<<<blocks, threads_per_block>>>(
          old_state->grid, old_state->x0_map, old_state->y0_map,
          old_state->radius_map, old_state->time_map);
      cudaError_t status = cudaGetLastError();
      if (status != cudaSuccess) {
        error = cuda_error("coordinate map kernel", status);
        return false;
      }
      return true;
    }
    release(*old_state);
    delete old_state;
    impl_ = nullptr;
  }

  int devices = 0;
  cudaError_t status = cudaGetDeviceCount(&devices);
  if (status != cudaSuccess)
    { error = cuda_error("cudaGetDeviceCount", status); return false; }
  if (devices <= 0) {
    error = "no CUDA device available";
    return false;
  }

  const std::size_t ncenters = static_cast<std::size_t>(grid.nx) * grid.ny;
  const std::size_t entries = ncenters * grid.nr * grid.nt *
                              grid.ne * grid.nphi;
  if (entries == 0 || entries > max_accumulator_entries) {
    error = "GPU Hough accumulator is too large";
    return false;
  }

  auto *state = new device_state_t;
  state->grid = grid;
  state->accumulator_entries = entries;
  status = cudaMalloc(&state->x0_map, entries * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc x0 map", status);
    release(*state); delete state; return false;
  }
  status = cudaMalloc(&state->y0_map, entries * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc y0 map", status);
    release(*state); delete state; return false;
  }
  status = cudaMalloc(&state->radius_map, entries * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc radius map", status);
    release(*state); delete state; return false;
  }
  status = cudaMalloc(&state->time_map, entries * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc time map", status);
    release(*state); delete state; return false;
  }
  status = cudaMalloc(&state->accumulator, entries * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc accumulator", status);
    release(*state); delete state;
    return false;
  }
  const std::size_t reduction_entries =
      (entries + threads_per_block - 1) / threads_per_block;
  status = cudaMalloc(&state->reduction_a,
                     reduction_entries * sizeof(unsigned long long));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc Hough reduction buffer", status);
    release(*state); delete state;
    return false;
  }
  state->reduction_capacity = reduction_entries;
  const int map_blocks = static_cast<int>(
      (entries + threads_per_block - 1) / threads_per_block);
  coordinate_map_kernel<<<map_blocks, threads_per_block>>>(
      state->grid, state->x0_map, state->y0_map,
      state->radius_map, state->time_map);
  status = cudaGetLastError();
  if (status != cudaSuccess) {
    error = cuda_error("coordinate map kernel", status);
    release(*state); delete state; return false;
  }
  impl_ = state;
  return true;
}

bool
engine_t::find(const float *x, const float *y, const float *time,
               int nhits, int max_candidates,
               candidate_t *candidates, int &candidate_count,
               std::string &error)
{
  auto *state = static_cast<device_state_t *>(impl_);
  if (!state) { error = "CUDA engine is not initialized"; return false; }
  if (nhits <= 0) { error = "cannot scan an empty event"; return false; }
  if (max_candidates <= 0) {
    error = "maximum candidate count must be positive";
    return false;
  }
  if (!candidates) {
    error = "candidate output buffer is null";
    return false;
  }
  candidate_count = 0;

  if (nhits > state->hit_capacity) {
    cudaFree(state->x); cudaFree(state->y); cudaFree(state->time);
    state->x = state->y = state->time = nullptr;
    cudaError_t status = cudaMalloc(&state->x, nhits * sizeof(float));
    if (status == cudaSuccess)
      status = cudaMalloc(&state->y, nhits * sizeof(float));
    if (status == cudaSuccess)
      status = cudaMalloc(&state->time, nhits * sizeof(float));
    if (status != cudaSuccess) {
      error = cuda_error("cudaMalloc event arrays", status);
      state->hit_capacity = 0;
      return false;
    }
    state->hit_capacity = nhits;
  }

  cudaError_t status = cudaMemcpy(state->x, x, nhits * sizeof(float),
                                  cudaMemcpyHostToDevice);
  if (status == cudaSuccess)
    status = cudaMemcpy(state->y, y, nhits * sizeof(float),
                        cudaMemcpyHostToDevice);
  if (status == cudaSuccess)
    status = cudaMemcpy(state->time, time, nhits * sizeof(float),
                        cudaMemcpyHostToDevice);
  if (status != cudaSuccess) {
    error = cuda_error("cudaMemcpy event arrays", status);
    return false;
  }

  const int blocks = static_cast<int>(
      (state->accumulator_entries + threads_per_block - 1) /
      threads_per_block);
  if (state->grid.circle_mode) {
    score_circle_kernel<<<blocks, threads_per_block>>>(
        state->x0_map, state->y0_map, state->radius_map, state->time_map,
        state->x, state->y, state->time, nhits, state->grid,
        state->accumulator);
  } else {
    score_kernel<<<blocks, threads_per_block>>>(
        state->x0_map, state->y0_map, state->radius_map, state->time_map,
        state->x, state->y, state->time, nhits, state->grid,
        state->accumulator);
  }
  status = cudaGetLastError();
  if (status == cudaSuccess)
    status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("Hough vote kernel", status);
    return false;
  }

  const std::size_t reduction_entries =
      (state->accumulator_entries + threads_per_block - 1) /
      threads_per_block;
  const int reduction_blocks = static_cast<int>(reduction_entries);
  reduce_block_max_kernel<<<reduction_blocks, threads_per_block>>>(
      state->accumulator, state->accumulator_entries, state->reduction_a);
  status = cudaGetLastError();
  if (status == cudaSuccess)
    status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("Hough maximum reduction", status);
    return false;
  }

  // One compact key is produced per GPU block. For the usual one-candidate
  // path, only this small array is copied back and the final comparison is
  // performed on the host, as in the Roberto implementation.
  std::vector<unsigned long long> block_maxima(reduction_entries);

  for (int candidate_number = 0;
       candidate_number < max_candidates; ++candidate_number) {
    status = cudaMemcpy(block_maxima.data(), state->reduction_a,
                        reduction_entries * sizeof(unsigned long long),
                        cudaMemcpyDeviceToHost);
    if (status != cudaSuccess) {
      error = cuda_error("cudaMemcpy Hough block maxima", status);
      return false;
    }
    const auto best = std::max_element(block_maxima.begin(),
                                       block_maxima.end());
    const unsigned long long best_key = *best;
    if (best_key == 0)
      break;
    const std::uint32_t score_bits =
        static_cast<std::uint32_t>(best_key >> 32);
    const std::uint32_t inverse_index =
        static_cast<std::uint32_t>(best_key & 0xffffffffull);
    const int best_index = static_cast<int>(0xffffffffu - inverse_index);
    float best_score = 0.f;
    std::memcpy(&best_score, &score_bits, sizeof(best_score));

    if (candidate_count == max_candidates)
      break;

    const int plane = state->grid.nr * state->grid.nt *
                      state->grid.ne * state->grid.nphi;
    const int center = best_index / plane;
    const int remainder = best_index % plane;
    const int ellipse_plane = state->grid.ne * state->grid.nphi;
    const int radius_plane = state->grid.nt * ellipse_plane;
    candidate_t candidate;
    candidate.x0_bin = center / state->grid.ny;
    candidate.y0_bin = center % state->grid.ny;
    candidate.radius_bin = remainder / radius_plane;
    const int radius_remainder = remainder % radius_plane;
    candidate.t_bin = radius_remainder / ellipse_plane;
    const int ellipse_remainder = radius_remainder % ellipse_plane;
    candidate.e_bin = ellipse_remainder / state->grid.nphi;
    candidate.phi_bin = ellipse_remainder % state->grid.nphi;
    candidate.score = best_score;
    candidates[candidate_count++] = candidate;

    if (candidate_number + 1 < max_candidates) {
      suppress_neighborhood_kernel<<<1, 1>>>(state->accumulator, state->grid,
                                             best_index);
      status = cudaGetLastError();
      if (status == cudaSuccess)
        status = cudaDeviceSynchronize();
      if (status != cudaSuccess) {
        error = cuda_error("Hough peak suppression", status);
        return false;
      }

      reduce_block_max_kernel<<<reduction_blocks, threads_per_block>>>(
          state->accumulator, state->accumulator_entries,
          state->reduction_a);
      status = cudaGetLastError();
      if (status == cudaSuccess)
        status = cudaDeviceSynchronize();
      if (status != cudaSuccess) {
        error = cuda_error("Hough block maximum reduction", status);
        return false;
      }
    }
  }
  return true;
}

ransac_engine_t::ransac_engine_t() = default;

ransac_engine_t::~ransac_engine_t()
{
  auto *state = static_cast<ransac_device_state_t *>(impl_);
  if (state) {
    release(*state);
    delete state;
  }
}

bool
ransac_engine_t::evaluate(
    const float *x, const float *y, const float *time, int nhits,
    ransac_model_t *models, int nmodels, double spatial_cut,
    double time_cut, bool write_masks, ransac_stats_t *stats,
    std::uint32_t *masks, std::string &error)
{
  if (nhits <= 0 || nmodels <= 0) {
    error = "RANSAC evaluation requires hits and hypotheses";
    return false;
  }
  if (!x || !y || !time || !models || !stats) {
    error = "RANSAC evaluation received a null buffer";
    return false;
  }
  if (write_masks && !masks) {
    error = "RANSAC inlier-mask buffer is null";
    return false;
  }

  auto *state = static_cast<ransac_device_state_t *>(impl_);
  if (!state) {
    int devices = 0;
    cudaError_t status = cudaGetDeviceCount(&devices);
    if (status != cudaSuccess) {
      error = cuda_error("cudaGetDeviceCount", status);
      return false;
    }
    if (devices <= 0) {
      error = "no CUDA device available";
      return false;
    }
    state = new ransac_device_state_t;
    impl_ = state;
  }

  cudaError_t status = cudaSuccess;
  if (nhits > state->hit_capacity) {
    cudaFree(state->x);
    cudaFree(state->y);
    cudaFree(state->time);
    state->x = state->y = state->time = nullptr;
    status = cudaMalloc(&state->x, nhits * sizeof(float));
    if (status == cudaSuccess)
      status = cudaMalloc(&state->y, nhits * sizeof(float));
    if (status == cudaSuccess)
      status = cudaMalloc(&state->time, nhits * sizeof(float));
    if (status != cudaSuccess) {
      error = cuda_error("cudaMalloc RANSAC event arrays", status);
      state->hit_capacity = 0;
      return false;
    }
    state->hit_capacity = nhits;
  }

  if (nmodels > state->model_capacity) {
    cudaFree(state->models);
    cudaFree(state->stats);
    state->models = nullptr;
    state->stats = nullptr;
    status = cudaMalloc(&state->models,
                        nmodels * sizeof(ransac_model_t));
    if (status == cudaSuccess)
      status = cudaMalloc(&state->stats,
                          nmodels * sizeof(ransac_stats_t));
    if (status != cudaSuccess) {
      error = cuda_error("cudaMalloc RANSAC hypothesis arrays", status);
      state->model_capacity = 0;
      return false;
    }
    state->model_capacity = nmodels;
  }

  const int mask_words = (nhits + 31) / 32;
  const std::size_t mask_entries =
      static_cast<std::size_t>(nmodels) * mask_words;
  if (write_masks && mask_entries > state->mask_capacity) {
    cudaFree(state->masks);
    state->masks = nullptr;
    status = cudaMalloc(&state->masks,
                        mask_entries * sizeof(std::uint32_t));
    if (status != cudaSuccess) {
      error = cuda_error("cudaMalloc RANSAC inlier masks", status);
      state->mask_capacity = 0;
      return false;
    }
    state->mask_capacity = mask_entries;
  }

  status = cudaMemcpy(state->x, x, nhits * sizeof(float),
                      cudaMemcpyHostToDevice);
  if (status == cudaSuccess)
    status = cudaMemcpy(state->y, y, nhits * sizeof(float),
                        cudaMemcpyHostToDevice);
  if (status == cudaSuccess)
    status = cudaMemcpy(state->time, time, nhits * sizeof(float),
                        cudaMemcpyHostToDevice);
  if (status == cudaSuccess)
    status = cudaMemcpy(state->models, models,
                        nmodels * sizeof(ransac_model_t),
                        cudaMemcpyHostToDevice);
  if (status != cudaSuccess) {
    error = cuda_error("cudaMemcpy RANSAC input", status);
    return false;
  }
  if (write_masks) {
    status = cudaMemset(state->masks, 0,
                        mask_entries * sizeof(std::uint32_t));
    if (status != cudaSuccess) {
      error = cuda_error("cudaMemset RANSAC inlier masks", status);
      return false;
    }
  }

  refine_ransac_kernel<<<nmodels, threads_per_block>>>(
      state->x, state->y, state->time, nhits, state->models, nmodels,
      spatial_cut, time_cut, write_masks, mask_words, state->stats,
      state->masks);
  status = cudaGetLastError();
  if (status == cudaSuccess)
    status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("RANSAC evaluation kernel", status);
    return false;
  }

  status = cudaMemcpy(stats, state->stats,
                      nmodels * sizeof(ransac_stats_t),
                      cudaMemcpyDeviceToHost);
  if (status == cudaSuccess)
    status = cudaMemcpy(models, state->models,
                        nmodels * sizeof(ransac_model_t),
                        cudaMemcpyDeviceToHost);
  if (status == cudaSuccess && write_masks)
    status = cudaMemcpy(masks, state->masks,
                        mask_entries * sizeof(std::uint32_t),
                        cudaMemcpyDeviceToHost);
  if (status != cudaSuccess) {
    error = cuda_error("cudaMemcpy RANSAC results", status);
    return false;
  }
  return true;
}

}

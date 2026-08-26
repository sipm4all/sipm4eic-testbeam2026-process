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
  unsigned long long *reduction_b = nullptr;
  int hit_capacity = 0;
  std::size_t accumulator_entries = 0;
  std::size_t reduction_capacity = 0;
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
reduce_scores_kernel(const float *scores, std::size_t entries,
                     unsigned long long *reduced)
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
    reduced[blockIdx.x] = keys[0];
}

__global__ void
reduce_keys_kernel(const unsigned long long *keys, std::size_t entries,
                   unsigned long long *reduced)
{
  __shared__ unsigned long long block_keys[threads_per_block];
  const int thread = threadIdx.x;
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) *
                            blockDim.x + thread;
  block_keys[thread] = index < entries ? keys[index] : 0;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (thread < stride)
      block_keys[thread] = max(block_keys[thread], block_keys[thread + stride]);
    __syncthreads();
  }
  if (thread == 0)
    reduced[blockIdx.x] = block_keys[0];
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
  cudaFree(state.reduction_b);
  state.reduction_a = nullptr;
  state.reduction_b = nullptr;
  state.hit_capacity = 0;
  state.accumulator_entries = 0;
  state.reduction_capacity = 0;
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
         first.ne == second.ne && first.nphi == second.nphi;
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
  status = cudaMalloc(&state->reduction_b,
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
  score_kernel<<<blocks, threads_per_block>>>(
      state->x0_map, state->y0_map, state->radius_map, state->time_map,
      state->x, state->y, state->time, nhits, state->grid,
      state->accumulator);
  status = cudaGetLastError();
  if (status == cudaSuccess)
    status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("Hough vote kernel", status);
    return false;
  }

  std::size_t reduction_entries =
      (state->accumulator_entries + threads_per_block - 1) /
      threads_per_block;
  int reduction_blocks = static_cast<int>(reduction_entries);
  reduce_scores_kernel<<<reduction_blocks, threads_per_block>>>(
      state->accumulator, state->accumulator_entries, state->reduction_a);
  status = cudaGetLastError();
  if (status == cudaSuccess)
    status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("Hough maximum reduction", status);
    return false;
  }

  unsigned long long *current = state->reduction_a;
  unsigned long long *next = state->reduction_b;
  for (std::size_t count = reduction_entries; count > 1;) {
    const std::size_t next_count = (count + threads_per_block - 1) /
                                   threads_per_block;
    reduce_keys_kernel<<<static_cast<int>(next_count), threads_per_block>>>(
        current, count, next);
    status = cudaGetLastError();
    if (status == cudaSuccess)
      status = cudaDeviceSynchronize();
    if (status != cudaSuccess) {
      error = cuda_error("Hough maximum reduction", status);
      return false;
    }
    std::swap(current, next);
    count = next_count;
  }

  for (int candidate_number = 0;
       candidate_number < max_candidates; ++candidate_number) {
    unsigned long long best_key = 0;
    status = cudaMemcpy(&best_key, current, sizeof(best_key),
                        cudaMemcpyDeviceToHost);
    if (status != cudaSuccess) {
      error = cuda_error("cudaMemcpy Hough maximum", status);
      return false;
    }
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

    suppress_neighborhood_kernel<<<1, 1>>>(state->accumulator, state->grid,
                                           best_index);
    status = cudaGetLastError();
    if (status == cudaSuccess)
      status = cudaDeviceSynchronize();
    if (status != cudaSuccess) {
      error = cuda_error("Hough peak suppression", status);
      return false;
    }

    reduce_scores_kernel<<<reduction_blocks, threads_per_block>>>(
        state->accumulator, state->accumulator_entries, state->reduction_a);
    status = cudaGetLastError();
    if (status == cudaSuccess)
      status = cudaDeviceSynchronize();
    if (status != cudaSuccess) {
      error = cuda_error("Hough maximum reduction", status);
      return false;
    }
    current = state->reduction_a;
    next = state->reduction_b;
    for (std::size_t count = reduction_entries; count > 1;) {
      const std::size_t next_count = (count + threads_per_block - 1) /
                                     threads_per_block;
      reduce_keys_kernel<<<static_cast<int>(next_count), threads_per_block>>>(
          current, count, next);
      status = cudaGetLastError();
      if (status == cudaSuccess)
        status = cudaDeviceSynchronize();
      if (status != cudaSuccess) {
        error = cuda_error("Hough maximum reduction", status);
        return false;
      }
      std::swap(current, next);
      count = next_count;
    }
  }
  return true;
}

}

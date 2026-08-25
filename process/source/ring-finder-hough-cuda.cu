#define RING_HOUGH_CUDA_IMPLEMENTATION
#include "ring-finder-hough-cuda.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
  int hit_capacity = 0;
  std::size_t accumulator_entries = 0;
  std::vector<float> host_accumulator;
};

struct ranked_entry_t {
  float score;
  int index;

  bool operator<(const ranked_entry_t &other) const
  {
    return score < other.score ||
           (score == other.score && index < other.index);
  }
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
                              grid.nr * grid.nt;
  if (index >= entries)
    return;

  const int plane = grid.nr * grid.nt;
  const int center = static_cast<int>(index / plane);
  const int remainder = static_cast<int>(index % plane);
  const int ix = center / grid.ny;
  const int iy = center % grid.ny;
  const int ir = remainder / grid.nt;
  const int it = remainder % grid.nt;

  x0_map[index] = static_cast<float>(grid.min_x0 + ix * grid.x0_step);
  y0_map[index] = static_cast<float>(grid.min_y0 + iy * grid.y0_step);
  radius_map[index] = static_cast<float>(grid.min_radius +
                                         ir * grid.radius_step);
  time_map[index] = static_cast<float>(grid.min_t + it * grid.t_step);
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
    const float dr = hypotf(x[hit] - x0, y[hit] - y0) - radius;
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
  state.hit_capacity = 0;
  state.accumulator_entries = 0;
  state.host_accumulator.clear();
}

bool
same_grid(const ring_hough_cuda::grid_t &first,
          const ring_hough_cuda::grid_t &second)
{
  return first.x0_step == second.x0_step &&
         first.y0_step == second.y0_step &&
         first.radius_step == second.radius_step &&
         first.t_step == second.t_step &&
         first.spatial_resolution == second.spatial_resolution &&
         first.time_resolution == second.time_resolution &&
         first.nx == second.nx && first.ny == second.ny &&
         first.nr == second.nr && first.nt == second.nt;
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
  const std::size_t entries = ncenters * grid.nr * grid.nt;
  if (entries == 0 || entries > max_accumulator_entries) {
    error = "GPU Hough accumulator is too large";
    return false;
  }

  auto *state = new device_state_t;
  state->grid = grid;
  state->accumulator_entries = entries;
  state->host_accumulator.resize(entries);
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
    delete state;
    return false;
  }
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

  status = cudaMemcpy(state->host_accumulator.data(), state->accumulator,
                      state->accumulator_entries * sizeof(float),
                      cudaMemcpyDeviceToHost);
  if (status != cudaSuccess) {
    error = cuda_error("cudaMemcpy Hough accumulator", status);
    return false;
  }

  std::vector<ranked_entry_t> ranked(state->host_accumulator.size());
  for (std::size_t index = 0; index < state->host_accumulator.size();
       ++index) {
    ranked[index].score = state->host_accumulator[index];
    ranked[index].index = static_cast<int>(index);
  }
  const std::size_t pool_limit = std::min(
      ranked.size(), static_cast<std::size_t>(max_candidates) * 64);
  if (pool_limit < ranked.size()) {
    std::nth_element(
        ranked.begin(), ranked.begin() + pool_limit, ranked.end(),
        [](const ranked_entry_t &first, const ranked_entry_t &second) {
          return second < first;
        });
    ranked.resize(pool_limit);
  }
  std::sort(ranked.begin(), ranked.end(),
            [](const ranked_entry_t &first, const ranked_entry_t &second) {
              return second < first;
            });

  // Find and suppress all requested maxima on the host after one accumulator
  // transfer. Repeating GPU reduction, host copies, and suppression kernels
  // for every candidate made peak extraction latency-dominated.
  for (const auto &entry : ranked) {
    if (candidate_count == max_candidates)
      break;
    const float best_score = state->host_accumulator[entry.index];
    if (best_score != entry.score)
      continue;
    if (!(best_score > 0.f))
      break;
    const int best_index = entry.index;

    const int plane = state->grid.nr * state->grid.nt;
    const int center = best_index / plane;
    const int remainder = best_index % plane;
    candidate_t candidate;
    candidate.x0_bin = center / state->grid.ny;
    candidate.y0_bin = center % state->grid.ny;
    candidate.radius_bin = remainder / state->grid.nt;
    candidate.t_bin = remainder % state->grid.nt;
    candidate.score = best_score;
    candidates[candidate_count++] = candidate;

    const int center_x = center / state->grid.ny;
    const int center_y = center % state->grid.ny;
    for (int dx = -1; dx <= 1; ++dx) {
      const int ix = center_x + dx;
      if (ix < 0 || ix >= state->grid.nx)
        continue;
      for (int dy = -1; dy <= 1; ++dy) {
        const int iy = center_y + dy;
        if (iy < 0 || iy >= state->grid.ny)
          continue;
        const int suppressed_center = ix * state->grid.ny + iy;
        for (int dr = -1; dr <= 1; ++dr) {
          const int ir = candidate.radius_bin + dr;
          if (ir < 0 || ir >= state->grid.nr)
            continue;
          for (int dt = -1; dt <= 1; ++dt) {
            const int it = candidate.t_bin + dt;
            if (it < 0 || it >= state->grid.nt)
              continue;
            const std::size_t index =
                static_cast<std::size_t>(suppressed_center) *
                    state->grid.nr * state->grid.nt +
                static_cast<std::size_t>(ir) * state->grid.nt + it;
            state->host_accumulator[index] = -1.f;
          }
        }
      }
    }
  }
  return true;
}

}

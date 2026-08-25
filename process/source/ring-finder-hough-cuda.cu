#define RING_HOUGH_CUDA_IMPLEMENTATION
#include "ring-finder-hough-cuda.h"

#include <cuda_runtime.h>

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
  float *x = nullptr;
  float *y = nullptr;
  float *time = nullptr;
  float *accumulator = nullptr;
  int *block_best_index = nullptr;
  float *block_best_score = nullptr;
  int hit_capacity = 0;
  int block_capacity = 0;
  int reduction_blocks = 0;
  std::size_t accumulator_entries = 0;
};

std::string
cuda_error(const char *operation, cudaError_t error)
{
  std::ostringstream message;
  message << operation << ": " << cudaGetErrorString(error);
  return message.str();
}

__global__ void
score_kernel(const float *x, const float *y, const float *time, int nhits,
             ring_hough_cuda::grid_t grid, float *accumulator)
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
  const float x0 = static_cast<float>(grid.min_x0 + ix * grid.x0_step);
  const float y0 = static_cast<float>(grid.min_y0 + iy * grid.y0_step);
  const float radius = static_cast<float>(grid.min_radius +
                                          ir * grid.radius_step);
  const float ring_time = static_cast<float>(grid.min_t + it * grid.t_step);
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

__global__ void
reduce_blocks_kernel(const float *accumulator, std::size_t entries,
                     int *block_best_index, float *block_best_score)
{
  extern __shared__ unsigned char storage[];
  auto *scores = reinterpret_cast<float *>(storage);
  auto *indices = reinterpret_cast<int *>(scores + blockDim.x);

  float local_score = -1.f;
  int local_index = -1;
  for (std::size_t index = static_cast<std::size_t>(blockIdx.x) *
                            blockDim.x + threadIdx.x;
       index < entries; index += static_cast<std::size_t>(gridDim.x) *
                                blockDim.x) {
    const float score = accumulator[index];
    if (score > local_score) {
      local_score = score;
      local_index = index;
    }
  }
  scores[threadIdx.x] = local_score;
  indices[threadIdx.x] = local_index;
  __syncthreads();

  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (threadIdx.x < stride && scores[threadIdx.x + stride] > scores[threadIdx.x]) {
      scores[threadIdx.x] = scores[threadIdx.x + stride];
      indices[threadIdx.x] = indices[threadIdx.x + stride];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    block_best_index[blockIdx.x] = indices[0];
    block_best_score[blockIdx.x] = scores[0];
  }
}

__global__ void
suppress_neighborhood_kernel(float *accumulator,
                             ring_hough_cuda::grid_t grid,
                             int center, int radius_bin, int time_bin)
{
  if (threadIdx.x != 0 || blockIdx.x != 0)
    return;
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
          const std::size_t index =
              static_cast<std::size_t>(suppressed_center) * grid.nr * grid.nt +
              static_cast<std::size_t>(ir) * grid.nt + it;
          accumulator[index] = -1.f;
        }
      }
    }
  }
}

void
release(device_state_t &state)
{
  cudaFree(state.x);
  cudaFree(state.y);
  cudaFree(state.time);
  cudaFree(state.accumulator);
  cudaFree(state.block_best_index);
  cudaFree(state.block_best_score);
  state = {};
}

bool
same_grid(const ring_hough_cuda::grid_t &first,
          const ring_hough_cuda::grid_t &second)
{
  return first.min_x0 == second.min_x0 &&
         first.min_y0 == second.min_y0 &&
         first.min_radius == second.min_radius &&
         first.min_t == second.min_t &&
         first.x0_step == second.x0_step &&
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
    if (same_grid(old_state->grid, grid))
      return true;
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
  status = cudaMalloc(&state->accumulator, entries * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc accumulator", status);
    delete state;
    return false;
  }
  state->reduction_blocks = static_cast<int>(
      (entries + threads_per_block - 1) / threads_per_block);
  state->block_capacity = state->reduction_blocks;
  status = cudaMalloc(&state->block_best_index,
                      state->block_capacity * sizeof(int));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc reduction indices", status);
    release(*state); delete state; return false;
  }
  status = cudaMalloc(&state->block_best_score,
                      state->block_capacity * sizeof(float));
  if (status != cudaSuccess) {
    error = cuda_error("cudaMalloc reduction scores", status);
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
      state->x, state->y, state->time, nhits, state->grid,
      state->accumulator);
  status = cudaGetLastError();
  if (status == cudaSuccess)
    status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("Hough vote kernel", status);
    return false;
  }

  for (int icandidate = 0; icandidate < max_candidates; ++icandidate) {
    reduce_blocks_kernel<<<state->reduction_blocks, threads_per_block,
                           threads_per_block * (sizeof(float) + sizeof(int))>>>(
        state->accumulator, state->accumulator_entries,
        state->block_best_index, state->block_best_score);
    status = cudaGetLastError();
    if (status == cudaSuccess)
      status = cudaDeviceSynchronize();
    if (status != cudaSuccess) {
      error = cuda_error("Hough reduction kernel", status);
      return false;
    }

    std::vector<int> block_indices(state->reduction_blocks);
    std::vector<float> block_scores(state->reduction_blocks);
    status = cudaMemcpy(block_indices.data(), state->block_best_index,
                        state->reduction_blocks * sizeof(int),
                        cudaMemcpyDeviceToHost);
    if (status == cudaSuccess)
      status = cudaMemcpy(block_scores.data(), state->block_best_score,
                          state->reduction_blocks * sizeof(float),
                          cudaMemcpyDeviceToHost);
    if (status != cudaSuccess) {
      error = cuda_error("cudaMemcpy Hough result", status);
      return false;
    }
    int best_index = -1;
    float best_score = 0.f;
    for (int iblock = 0; iblock < state->reduction_blocks; ++iblock) {
      if (block_scores[iblock] > best_score) {
        best_score = block_scores[iblock];
        best_index = block_indices[iblock];
      }
    }
    if (best_index < 0 || !(best_score > 0.f))
      break;

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

    suppress_neighborhood_kernel<<<1, 1>>>(
        state->accumulator, state->grid, center,
        candidate.radius_bin, candidate.t_bin);
    status = cudaGetLastError();
    if (status != cudaSuccess) {
      error = cuda_error("Hough peak suppression kernel", status);
      return false;
    }
  }
  status = cudaDeviceSynchronize();
  if (status != cudaSuccess) {
    error = cuda_error("Hough peak suppression", status);
    return false;
  }
  return true;
}

}

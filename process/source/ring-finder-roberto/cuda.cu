#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "common.h"

namespace {

constexpr int threads_per_block = 256;
constexpr int point_capacity = 1 << 16;

void
handle_error(cudaError_t error, const char *file, int line)
{
  if (error != cudaSuccess) {
    std::fprintf(stderr, "%s in %s at line %d\n",
                 cudaGetErrorString(error), file, line);
    std::exit(EXIT_FAILURE);
  }
}

#define HANDLE_ERROR(error) handle_error((error), __FILE__, __LINE__)

data_t d_data{};
int allocated_size = 0;

__global__ void
hough_init_kernel(data_t data)
{
  const int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= size)
    return;

  const int ix = index % data.bins.x;
  const int iy = (index / data.bins.x) % data.bins.y;
  const int ir = (index / (data.bins.x * data.bins.y)) % data.bins.r;
  const int it = (index / (data.bins.x * data.bins.y * data.bins.r)) %
                data.bins.t;

  data.map.x[index] = data.min.x +
                      (0.5f + ix) * (data.max.x - data.min.x) / data.bins.x;
  data.map.y[index] = data.min.y +
                      (0.5f + iy) * (data.max.y - data.min.y) / data.bins.y;
  data.map.r[index] = data.min.r +
                      (0.5f + ir) * (data.max.r - data.min.r) / data.bins.r;
  data.map.t[index] = data.min.t +
                      (0.5f + it) * (data.max.t - data.min.t) / data.bins.t;
}

__global__ void
hough_transform_kernel(data_t data)
{
  const int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= size)
    return;

  const float center_x = data.map.x[index];
  const float center_y = data.map.y[index];
  const float radius = data.map.r[index];
  const float ring_time = data.map.t[index];
  const float time_sigma = data.sigma.t;
  const float time_exponent = -0.5f / (time_sigma * time_sigma);

  data.hough.h[index] = 0.f;
  data.hough.nh[index] = 0;
  for (int hit = 0; hit < data.points.n; ++hit) {
    if (data.points.x[hit] == 0.f || data.points.y[hit] == 0.f)
      continue;

    const float dx = data.points.x[hit] - center_x;
    const float dy = data.points.y[hit] - center_y;
    const float dt = data.points.t[hit] - ring_time;
    const float dr = hypotf(dx, dy) - radius;
    const float spatial_weight = expf(-0.11337868f * dr * dr);
    const float time_weight = expf(time_exponent * dt * dt);
    data.hough.h[index] += spatial_weight * time_weight;
    ++data.hough.nh[index];
  }
}

__global__ void
find_max_kernel(data_t data)
{
  __shared__ float scores[threads_per_block];
  __shared__ int indices[threads_per_block];

  const int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  const int global_index = blockIdx.x * blockDim.x + threadIdx.x;
  const int thread = threadIdx.x;
  scores[thread] = global_index < size ? data.hough.h[global_index] : -INFINITY;
  indices[thread] = global_index < size ? global_index : -1;
  __syncthreads();

  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (thread < stride && scores[thread + stride] > scores[thread]) {
      scores[thread] = scores[thread + stride];
      indices[thread] = indices[thread + stride];
    }
    __syncthreads();
  }

  if (thread == 0) {
    data.hough.rh[blockIdx.x] = scores[0];
    data.hough.rhi[blockIdx.x] = indices[0];
  }
}

}

void
hough_init(data_t host_data)
{
  const int size = host_data.bins.x * host_data.bins.y *
                   host_data.bins.r * host_data.bins.t;
  const int grid_size = (size + threads_per_block - 1) / threads_per_block;

  if (allocated_size != 0 && allocated_size != size) {
    hough_free();
  }

  if (allocated_size == 0) {
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.map.x),
                            size * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.map.y),
                            size * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.map.r),
                            size * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.map.t),
                            size * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.points.x),
                            point_capacity * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.points.y),
                            point_capacity * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.points.t),
                            point_capacity * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.hough.h),
                            size * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.hough.nh),
                            size * sizeof(int)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.hough.rh),
                            grid_size * sizeof(float)));
    HANDLE_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_data.hough.rhi),
                            grid_size * sizeof(int)));
    allocated_size = size;
  }

  d_data.min = host_data.min;
  d_data.max = host_data.max;
  d_data.bins = host_data.bins;
  d_data.sigma = host_data.sigma;
  hough_init_kernel<<<grid_size, threads_per_block>>>(d_data);
  HANDLE_ERROR(cudaGetLastError());
  HANDLE_ERROR(cudaDeviceSynchronize());

  HANDLE_ERROR(cudaMemcpy(host_data.map.x, d_data.map.x,
                          size * sizeof(float), cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_data.map.y, d_data.map.y,
                          size * sizeof(float), cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_data.map.r, d_data.map.r,
                          size * sizeof(float), cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_data.map.t, d_data.map.t,
                          size * sizeof(float), cudaMemcpyDeviceToHost));
}

void
hough_free()
{
  HANDLE_ERROR(cudaFree(d_data.map.x));
  HANDLE_ERROR(cudaFree(d_data.map.y));
  HANDLE_ERROR(cudaFree(d_data.map.r));
  HANDLE_ERROR(cudaFree(d_data.map.t));
  HANDLE_ERROR(cudaFree(d_data.points.x));
  HANDLE_ERROR(cudaFree(d_data.points.y));
  HANDLE_ERROR(cudaFree(d_data.points.t));
  HANDLE_ERROR(cudaFree(d_data.hough.h));
  HANDLE_ERROR(cudaFree(d_data.hough.nh));
  HANDLE_ERROR(cudaFree(d_data.hough.rh));
  HANDLE_ERROR(cudaFree(d_data.hough.rhi));
  d_data = {};
  allocated_size = 0;
}

void
hough_transform(data_t host_data)
{
  const int size = host_data.bins.x * host_data.bins.y *
                   host_data.bins.r * host_data.bins.t;
  const int grid_size = (size + threads_per_block - 1) / threads_per_block;

  HANDLE_ERROR(cudaMemcpy(d_data.points.x, host_data.points.x,
                          host_data.points.n * sizeof(float),
                          cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(d_data.points.y, host_data.points.y,
                          host_data.points.n * sizeof(float),
                          cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(d_data.points.t, host_data.points.t,
                          host_data.points.n * sizeof(float),
                          cudaMemcpyHostToDevice));
  d_data.points.n = host_data.points.n;

  hough_transform_kernel<<<grid_size, threads_per_block>>>(d_data);
  HANDLE_ERROR(cudaGetLastError());
  find_max_kernel<<<grid_size, threads_per_block>>>(d_data);
  HANDLE_ERROR(cudaGetLastError());
  HANDLE_ERROR(cudaDeviceSynchronize());

  HANDLE_ERROR(cudaMemcpy(host_data.hough.h, d_data.hough.h,
                          size * sizeof(float), cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_data.hough.nh, d_data.hough.nh,
                          size * sizeof(int), cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_data.hough.rh, d_data.hough.rh,
                          grid_size * sizeof(float), cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_data.hough.rhi, d_data.hough.rhi,
                          grid_size * sizeof(int), cudaMemcpyDeviceToHost));
}

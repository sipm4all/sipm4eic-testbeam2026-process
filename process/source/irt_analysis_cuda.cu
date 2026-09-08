#include "irt_analysis_cuda.h"
#include <cuda_runtime.h>
#include <cmath>

namespace irt_analysis_cuda {
constexpr int NT=100, NX=64; constexpr float TMIN=0.f, TMAX=.1f, XMIN=-32.f, XMAX=32.f;
struct buffers_t { float *th=nullptr,*tm=nullptr,*score=nullptr; std::size_t capacity=0; };
static buffers_t buffers;
static void release() { cudaFree(buffers.th); cudaFree(buffers.tm); cudaFree(buffers.score); buffers={}; }
static bool reserve(std::size_t n) {
  if (n <= buffers.capacity) return true;
  release();
  if (cudaMalloc(&buffers.th,n*sizeof(float))!=cudaSuccess ||
      cudaMalloc(&buffers.tm,n*sizeof(float))!=cudaSuccess ||
      cudaMalloc(&buffers.score,NT*NX*sizeof(float))!=cudaSuccess) { release(); return false; }
  buffers.capacity=n; return true;
}
__global__ void score(const float *th, const float *tm, int n, float *out) {
  const int cell=blockIdx.x*blockDim.x+threadIdx.x; if(cell>=NT*NX)return;
  const int it=cell/NX, ix=cell%NX;
  const float t=TMIN+(it+.5f)*(TMAX-TMIN)/NT, x=XMIN+(ix+.5f)*(XMAX-XMIN)/NX;
  float s=0.f;
  for(int i=0;i<n;++i){float a=(th[i]-t)/.002f,b=(tm[i]-x);s+=expf(-.5f*(a*a+b*b));}
  out[cell]=s;
}
bool available(){int n=0; return cudaGetDeviceCount(&n)==cudaSuccess&&n>0;}
bool peak(const std::vector<float>& th,const std::vector<float>& tm,peak_t& r){
  if(th.empty()||th.size()!=tm.size() || !reserve(th.size()))return false;
  cudaMemcpy(buffers.th,th.data(),th.size()*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(buffers.tm,tm.data(),tm.size()*sizeof(float),cudaMemcpyHostToDevice);
  score<<<(NT*NX+255)/256,256>>>(buffers.th,buffers.tm,th.size(),buffers.score); cudaDeviceSynchronize();
  float values[NT*NX]; cudaMemcpy(values,buffers.score,sizeof(values),cudaMemcpyDeviceToHost);
  int best=0; for(int i=1;i<NT*NX;++i)if(values[i]>values[best])best=i;
  r.theta=TMIN+(best/NX+.5f)*(TMAX-TMIN)/NT; r.time=XMIN+(best%NX+.5f)*(XMAX-XMIN)/NX; r.score=values[best]; return true;
}
}

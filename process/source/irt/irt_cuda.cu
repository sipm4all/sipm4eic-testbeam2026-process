#include "irt_cuda.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>

namespace irt::cuda {
namespace {

struct v3 { double x, y, z; };

__device__ v3 add(v3 a, v3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
__device__ v3 sub(v3 a, v3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
__device__ v3 mul(v3 a, double s) { return {a.x*s, a.y*s, a.z*s}; }
__device__ double dot(v3 a, v3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
__device__ v3 cross(v3 a, v3 b) {
  return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
__device__ double norm(v3 a) { return sqrt(dot(a,a)); }
__device__ v3 unit(v3 a) { double n=norm(a); return mul(a, 1./n); }
__device__ v3 rotate(v3 v, double angle, v3 axis) {
  axis = unit(axis);
  return add(add(mul(v, cos(angle)), mul(cross(axis,v), sin(angle))),
             mul(axis, dot(axis,v)*(1.-cos(angle))));
}

__device__ bool inverse_ray(v3 emission, v3 detector, v3 center,
                            double radius, v3 &reflection) {
  const v3 ce=sub(emission,center), cd=sub(detector,center);
  const double a=norm(ce), d=norm(cd), alpha=acos(dot(ce,cd)/(a*d));
  const double sine=sin(alpha);
  if (!(a>0.) || !(d>0.) || !(radius>0.) || fabs(sine)<1.e-10) return false;
  auto f = [&](double beta) {
    return a*d*sin(alpha-2.*beta)+radius*(a*sin(beta)-d*sin(alpha-beta));
  };
  auto df = [&](double beta) {
    return -2.*a*d*cos(alpha-2.*beta)+radius*(a*cos(beta)+d*cos(alpha-beta));
  };
  double beta=.5*alpha;
  for (int i=0;i<100;++i) {
    const double der=df(beta); if (!isfinite(der) || fabs(der)<1.e-12) return false;
    const double delta=f(beta)/der; if (!isfinite(delta) || fabs(delta)>1.) return false;
    beta-=delta; if (fabs(delta)<1.e-10) break;
  }
  reflection=add(center, add(mul(ce, radius*cos(beta)/a-radius*sin(beta)*cos(alpha)/(a*sine)),
                         mul(cd, radius*sin(beta)/(d*sine))));
  return isfinite(reflection.x)&&isfinite(reflection.y)&&isfinite(reflection.z);
}

__device__ bool reconstruct_one(const geometry_t &g, float xin, float yin,
                                float tin, float &theta, float &phi,
                                float &etime) {
  v3 normal=unit({g.detector_normal[0],g.detector_normal[1],g.detector_normal[2]});
  v3 ref={0.,1.,0.}; if (fabs(dot(ref,normal))>.999) ref={1.,0.,0.};
  v3 ax=unit(sub(ref,mul(normal,dot(ref,normal))));
  v3 ay=cross(normal,ax);
  const v3 plane_normal = normal;
  v3 detector_rotation={g.detector_rotation[0],g.detector_rotation[1],g.detector_rotation[2]};
  if (norm(detector_rotation)>0.)
    normal=rotate(normal,norm(detector_rotation),detector_rotation);
  normal=rotate(normal,g.detector_tilt_x,ax);
  normal=rotate(normal,g.detector_tilt_y,ay);
  v3 plane=add({g.detector_center[0],g.detector_center[1],g.detector_center[2]},
               mul(plane_normal,g.detector_radius));
  v3 local={xin,yin,0.};
  // RotateUz equivalent for a vector initially in the detector xy plane.
  // Exact equivalent of ROOT's TVector3::RotateUz for local z = 0.
  const double up = hypot(normal.x, normal.y);
  if (up > 1.e-15) {
    local = {(normal.x * normal.z * xin - normal.y * yin) / up,
             (normal.y * normal.z * xin + normal.x * yin) / up,
             -up * xin};
  } else if (normal.z < 0.) {
    local = {-xin, yin, 0.};
  }
  plane=add(plane,local);
  v3 mc={g.mirror_center[0],g.mirror_center[1],g.mirror_center[2]};
  v3 mr={g.mirror_rotation[0],g.mirror_rotation[1],g.mirror_rotation[2]};
  if (norm(mr)>0.) mc=add(g.mirror_pivot[0] == g.mirror_pivot[0] ?
      rotate(sub(mc,{g.mirror_pivot[0],g.mirror_pivot[1],g.mirror_pivot[2]}),norm(mr),mr) : mc,
      {g.mirror_pivot[0],g.mirror_pivot[1],g.mirror_pivot[2]});
  v3 reflection;
  if (!inverse_ray({g.emission[0],g.emission[1],g.emission[2]},plane,mc,g.mirror_radius,reflection)) return false;
  v3 direction=unit(sub(reflection,{g.emission[0],g.emission[1],g.emission[2]}));
  v3 track={g.track[0],g.track[1],g.track[2]};
  if (!isfinite(norm(direction)) || dot(direction,track)<=0.) return false;
  theta=static_cast<float>(acos(fmax(-1.,fmin(1.,dot(track,direction)))));
  v3 tr={0.,1.,0.}; if (fabs(dot(tr,track))>.999) tr={1.,0.,0.};
  v3 tx=unit(sub(tr,mul(track,dot(tr,track))));
  v3 ty=cross(track,tx);
  phi=static_cast<float>(atan2(dot(direction,ty),dot(direction,tx)));
  const double path=norm(sub(reflection,plane))+norm(sub(reflection,{g.emission[0],g.emission[1],g.emission[2]}));
  etime=static_cast<float>(tin*3.125-path/299.792458);
  return isfinite(theta)&&isfinite(phi)&&isfinite(etime);
}

__global__ void reconstruct_kernel(geometry_t g, const float *x, const float *y,
                                    const float *time, int n, float *theta,
                                    float *phi, float *etime, unsigned char *valid) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x; if (i>=n) return;
  valid[i]=reconstruct_one(g,x[i],y[i],time[i],theta[i],phi[i],etime[i]) ? 1 : 0;
  if (!valid[i]) { theta[i]=NAN; phi[i]=NAN; etime[i]=NAN; }
}

} // namespace

bool available() { int n=0; return cudaGetDeviceCount(&n)==cudaSuccess && n>0; }

bool reconstruct(const geometry_t &g, const float *x, const float *y,
                 const float *time, int n, float *theta, float *phi,
                 float *etime, unsigned char *valid) {
  if (n <= 0) return true;
  if (!available()) return false;
  float *dx=nullptr,*dy=nullptr,*dt=nullptr,*do1=nullptr,*dp=nullptr,*de=nullptr;
  unsigned char *dv=nullptr; const size_t bytes=static_cast<size_t>(n)*sizeof(float);
  auto ok = [&](cudaError_t e, const char *where){
    if (e == cudaSuccess) return true;
    fprintf(stderr, "CUDA IRT error at %s: %s\\n", where, cudaGetErrorString(e));
    return false;
  };
  if (!ok(cudaMalloc(&dx,bytes),"malloc x")||!ok(cudaMalloc(&dy,bytes),"malloc y")||
      !ok(cudaMalloc(&dt,bytes),"malloc time")||!ok(cudaMalloc(&do1,bytes),"malloc theta")||
      !ok(cudaMalloc(&dp,bytes),"malloc phi")||!ok(cudaMalloc(&de,bytes),"malloc emission time")||
      !ok(cudaMalloc(&dv,static_cast<size_t>(n)),"malloc valid")) return false;
  bool result=ok(cudaMemcpy(dx,x,bytes,cudaMemcpyHostToDevice),"copy x")&&
    ok(cudaMemcpy(dy,y,bytes,cudaMemcpyHostToDevice),"copy y")&&
    ok(cudaMemcpy(dt,time,bytes,cudaMemcpyHostToDevice),"copy time");
  if (result) {
    reconstruct_kernel<<<(n+255)/256,256>>>(g,dx,dy,dt,n,do1,dp,de,dv);
    result=ok(cudaGetLastError(),"launch")&&ok(cudaDeviceSynchronize(),"synchronize");
  }
  if (result) result=ok(cudaMemcpy(theta,do1,bytes,cudaMemcpyDeviceToHost),"copy theta")&&
    ok(cudaMemcpy(phi,dp,bytes,cudaMemcpyDeviceToHost),"copy phi")&&
    ok(cudaMemcpy(etime,de,bytes,cudaMemcpyDeviceToHost),"copy emission time")&&
    ok(cudaMemcpy(valid,dv,static_cast<size_t>(n),cudaMemcpyDeviceToHost),"copy valid");
  cudaFree(dx);cudaFree(dy);cudaFree(dt);cudaFree(do1);cudaFree(dp);cudaFree(de);cudaFree(dv); return result;
}

} // namespace irt::cuda

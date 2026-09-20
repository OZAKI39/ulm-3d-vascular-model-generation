#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>
#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){std::fprintf(stderr,"%s: %s\n",#call,cudaGetErrorString(e));return 1;} } while(0)
__global__ void twice_plus_one(const int *a,int *b,int n) {
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<n)b[i]=2*a[i]+1;
}
int main() {
  int count=0,runtime=0,driver=0;cudaDeviceProp p;
  CHECK(cudaGetDeviceCount(&count));if(count<1)return 2;
  CHECK(cudaSetDevice(0));CHECK(cudaGetDeviceProperties(&p,0));
  CHECK(cudaRuntimeGetVersion(&runtime));CHECK(cudaDriverGetVersion(&driver));
  if(!std::strstr(p.name,"RTX 4090") || p.major!=8 || p.minor!=9)return 3;
  if(__CUDACC_VER_MAJOR__!=12 || __CUDACC_VER_MINOR__!=6 || runtime!=12060)return 4;
  const int n=1024;int input[n],output[n],*device_in,*device_out;
  for(int i=0;i<n;++i)input[i]=i;
  CHECK(cudaMalloc(&device_in,sizeof(input)));CHECK(cudaMalloc(&device_out,sizeof(output)));
  CHECK(cudaMemcpy(device_in,input,sizeof(input),cudaMemcpyHostToDevice));
  twice_plus_one<<<8,128>>>(device_in,device_out,n);
  CHECK(cudaGetLastError());CHECK(cudaDeviceSynchronize());
  CHECK(cudaMemcpy(output,device_out,sizeof(output),cudaMemcpyDeviceToHost));
  for(int i=0;i<n;++i)if(output[i]!=2*i+1)return 5;
  CHECK(cudaFree(device_in));CHECK(cudaFree(device_out));CHECK(cudaGetLastError());
  std::printf("CUDA_SMOKE count=%d GPU=%s arch=%d.%d runtime=%d driver_api=%d nvcc=%d.%d.%d elements=%d correct=1 last_error=0\n",count,p.name,p.major,p.minor,runtime,driver,__CUDACC_VER_MAJOR__,__CUDACC_VER_MINOR__,__CUDACC_VER_BUILD__,n);
  return 0;
}

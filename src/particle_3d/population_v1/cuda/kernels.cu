#include <torch/extension.h>
#include <ATen/cuda/CUDAContext.h>
#include <c10/cuda/CUDAGuard.h>
#include <c10/cuda/CUDAException.h>
#include <cuda_runtime.h>
#include <vector>
#include <cmath>
#include <climits>
#include <cfloat>

using i64 = int64_t;

__device__ void weights(const double* p, i64 cell, const double* origin,
                       const double* inverse, double* w) {
  double dx[3];for(int k=0;k<3;k++)dx[k]=p[k]-origin[3*cell+k];
  const double* m=inverse+9*cell;
  // numpy einsum's short 3-vector accumulation uses (term0+term2)+term1.
  for(int k=0;k<3;k++)w[k+1]=(m[3*k]*dx[0]+m[3*k+2]*dx[2])+m[3*k+1]*dx[1];
  w[0]=1.-((w[1]+w[2])+w[3]);
}
__device__ bool contains(const double* w,double tol) {
  double low=w[0],high=w[0];for(int k=1;k<4;k++){low=fmin(low,w[k]);high=fmax(high,w[k]);}
  return low>=-tol && high<=1.+tol;
}
__device__ bool box_contains(const double* p,const double* box,double pad) {
  for(int k=0;k<3;k++)if(p[k]<box[k]-pad||p[k]>box[k+3]+pad)return false;
  return true;
}
__device__ i64 bvh_locate(const double* p,const double* box,const i64* link,
    const i64* order,const double* origin,const double* inverse,const double* tau,double pad,double* result) {
  int stack[128];int top=0;stack[top++]=0;i64 best=LLONG_MAX;
  while(top){int node=stack[--top];if(!box_contains(p,box+6*node,pad))continue;
    const i64* l=link+4*node;
    if(l[0]>=0){stack[top++]=l[1];stack[top++]=l[0];}
    else for(i64 z=l[2];z<l[2]+l[3];z++){
      i64 c=order[z];if(c>=best)continue;double w[4];weights(p,c,origin,inverse,w);
      if(contains(w,tau[c])){best=c;for(int k=0;k<4;k++)result[k]=w[k];}
    }
  }
  return best==LLONG_MAX?-1:best;
}
__global__ void fem_kernel(const double* x,const i64* hints,i64 count,i64 cell_count,
    const i64* tetra,const double* origin,const double* inverse,const double* tau,const i64* adjacent,
    const double* shape_norm,const double* box,const i64* link,const i64* order,
    const double* velocity,const double* pressure,const double* gradient,const double* curl,
    const double* strain,double pad,i64* ids,double* output,i64* stats) {
  i64 row=blockIdx.x*blockDim.x+threadIdx.x;if(row>=count)return;
  const double* p=x+3*row;double w[4];i64 cell=hints[row];bool local=false;int walks=0;
  for(int k=0;k<64 && cell>=0 && cell<cell_count;k++){
    weights(p,cell,origin,inverse,w);
    if(contains(w,tau[cell])){
      bool unique=true;
      for(int j=0;j<4;j++)if(w[j]/shape_norm[4*cell+j]<=pad)unique=false;
      // At a face/edge/vertex all tolerance-containing cells must be considered.
      if(unique)local=true;
      break;
    }
    int face=0;for(int j=1;j<4;j++)if(w[j]<w[face])face=j;
    cell=adjacent[4*cell+face];walks++;
  }
  if(!local)cell=bvh_locate(p,box,link,order,origin,inverse,tau,pad,w);
  ids[row]=cell;stats[2*row]=walks;stats[2*row+1]=!local;
  double* out=output+25*row;
  if(cell<0){for(int k=0;k<25;k++)out[k]=NAN;return;}
  const i64* nodes=tetra+4*cell;
  // Fixed-order interpolation, no atomics and no reduced precision.
  for(int k=0;k<3;k++){
    out[k]=((w[0]*velocity[3*nodes[0]+k]+w[1]*velocity[3*nodes[1]+k])+
             w[2]*velocity[3*nodes[2]+k])+w[3]*velocity[3*nodes[3]+k];
  }
  out[3]=((w[0]*pressure[nodes[0]]+w[1]*pressure[nodes[1]])+w[2]*pressure[nodes[2]])+w[3]*pressure[nodes[3]];
  for(int k=0;k<9;k++){out[4+k]=gradient[9*cell+k];out[16+k]=strain[9*cell+k];}
  for(int k=0;k<3;k++)out[13+k]=curl[3*cell+k];
}
std::vector<torch::Tensor> fem_sample_cuda(std::vector<torch::Tensor> a,double padding) {
  TORCH_CHECK(a.size()==16,"FEM tensor bundle size");
  for(auto& t:a)TORCH_CHECK(t.is_cuda()&&t.is_contiguous(),"CUDA contiguous inputs required");
  TORCH_CHECK(a[0].scalar_type()==torch::kFloat64 && a[1].scalar_type()==torch::kInt64,"FP64 coordinates/int64 IDs required");
  c10::cuda::CUDAGuard guard(a[0].device());i64 n=a[0].size(0);
  auto ids=torch::empty({n},a[1].options());auto result=torch::empty({n,25},a[0].options());
  auto stats=torch::empty({n,2},a[1].options());
  if(n)fem_kernel<<<(n+127)/128,128,0,at::cuda::getCurrentCUDAStream()>>>(
    a[0].data_ptr<double>(),a[1].data_ptr<i64>(),n,a[2].size(0),
    a[2].data_ptr<i64>(),a[3].data_ptr<double>(),a[4].data_ptr<double>(),a[5].data_ptr<double>(),
    a[6].data_ptr<i64>(),a[7].data_ptr<double>(),a[8].data_ptr<double>(),a[9].data_ptr<i64>(),a[10].data_ptr<i64>(),
    a[11].data_ptr<double>(),a[12].data_ptr<double>(),a[13].data_ptr<double>(),a[14].data_ptr<double>(),a[15].data_ptr<double>(),
    padding,ids.data_ptr<i64>(),result.data_ptr<double>(),stats.data_ptr<i64>());
  C10_CUDA_KERNEL_LAUNCH_CHECK();return {ids,result,stats};
}

__device__ i64 cell_key(i64 x,i64 y,i64 z,i64 real) {
  uint64_t key=(uint64_t(x)*73856093ULL)^(uint64_t(y)*19349663ULL)^
    (uint64_t(z)*83492791ULL)^(uint64_t(real)*2654435761ULL);
  return i64(key&0x7fffffffffffffffULL);
}
__device__ i64 lower_key(const i64* keys,i64 n,i64 key) {
  i64 l=0,r=n;while(l<r){i64 m=(l+r)/2;if(keys[m]<key)l=m+1;else r=m;}return l;
}
template<bool fill> __global__ void neighbor_kernel(const double* x,const double* radii,
    const double* bound,const i64* realization,const i64* cells,const i64* keys,
    const i64* sorted_index,i64 n,double size,i64* counts,const i64* offsets,
    i64* pairs,double* geometry) {
  i64 i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;i64 count=0;
  i64 cx=cells[3*i],cy=cells[3*i+1],cz=cells[3*i+2],rid=realization[i];
  for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
    i64 ax=cx+dx,ay=cy+dy,az=cz+dz,key=cell_key(ax,ay,az,rid);
    for(i64 k=lower_key(keys,n,key);k<n && keys[k]==key;k++){
      i64 j=sorted_index[k];
      if(j<=i || realization[j]!=rid || cells[3*j]!=ax || cells[3*j+1]!=ay || cells[3*j+2]!=az)continue;
      double rx=x[3*i]-x[3*j],ry=x[3*i+1]-x[3*j+1],rz=x[3*i+2]-x[3*j+2];
      double distance=sqrt((rx*rx+ry*ry)+rz*rz);
      if(distance>bound[i]+bound[j])continue;
      if(fill){i64 z=offsets[i]+count;pairs[2*z]=i;pairs[2*z+1]=j;
        geometry[4*z]=distance-radii[i]-radii[j];
        geometry[4*z+1]=distance?rx/distance:1.;geometry[4*z+2]=distance?ry/distance:0.;geometry[4*z+3]=distance?rz/distance:0.;
      }count++;
    }
  }
  if(!fill)counts[i]=count;
}
torch::Tensor neighbor_count_cuda(std::vector<torch::Tensor> a,double size) {
  TORCH_CHECK(a.size()==7,"Neighbor bundle size");c10::cuda::CUDAGuard guard(a[0].device());
  for(auto& t:a)TORCH_CHECK(t.is_cuda()&&t.is_contiguous(),"Contiguous device arrays required");
  i64 n=a[0].size(0);auto counts=torch::zeros({n},a[3].options());
  if(n)neighbor_kernel<false><<<(n+127)/128,128,0,at::cuda::getCurrentCUDAStream()>>>(
    a[0].data_ptr<double>(),a[1].data_ptr<double>(),a[2].data_ptr<double>(),a[3].data_ptr<i64>(),
    a[4].data_ptr<i64>(),a[5].data_ptr<i64>(),a[6].data_ptr<i64>(),n,size,counts.data_ptr<i64>(),nullptr,nullptr,nullptr);
  C10_CUDA_KERNEL_LAUNCH_CHECK();return counts;
}
std::vector<torch::Tensor> neighbor_fill_cuda(std::vector<torch::Tensor> a,double size,torch::Tensor offsets,i64 total) {
  c10::cuda::CUDAGuard guard(a[0].device());i64 n=a[0].size(0);
  auto pairs=torch::empty({total,2},a[3].options());auto geometry=torch::empty({total,4},a[0].options());
  if(n)neighbor_kernel<true><<<(n+127)/128,128,0,at::cuda::getCurrentCUDAStream()>>>(
    a[0].data_ptr<double>(),a[1].data_ptr<double>(),a[2].data_ptr<double>(),a[3].data_ptr<i64>(),
    a[4].data_ptr<i64>(),a[5].data_ptr<i64>(),a[6].data_ptr<i64>(),n,size,nullptr,offsets.data_ptr<i64>(),pairs.data_ptr<i64>(),geometry.data_ptr<double>());
  C10_CUDA_KERNEL_LAUNCH_CHECK();return {pairs,geometry};
}

__device__ double dot3(const double* a,const double* b){return (a[0]*b[0]+a[1]*b[1])+a[2]*b[2];}
__device__ double ein3(const double* a,const double* b){return (a[0]*b[0]+a[2]*b[2])+a[1]*b[1];}
__device__ void triangle_closest(const double* p,const double* a,double* best,double* bw) {
  double e[3],f[3],v[3],d[3];for(int j=0;j<3;j++){e[j]=a[3+j]-a[j];f[j]=a[6+j]-a[j];v[j]=p[j]-a[j];}
  double ee=ein3(e,e),ef=ein3(e,f),ff=ein3(f,f),ev=ein3(e,v),fv=ein3(f,v),den=ee*ff-ef*ef;
  double u=den>0?(ff*ev-ef*fv)/den:0.,w=den>0?(ee*fv-ef*ev)/den:0.;
  bw[0]=1.-u-w;bw[1]=u;bw[2]=w;
  for(int j=0;j<3;j++){best[j]=(a[j]+u*e[j])+w*f[j];d[j]=p[j]-best[j];}
  double distance=(den>0&&bw[0]>=0&&u>=0&&w>=0)?dot3(d,d):INFINITY;
  for(int edge=0;edge<3;edge++){
    int i=edge,j=(edge+1)%3;double ee2[3],vv[3],cp[3];
    for(int h=0;h<3;h++){ee2[h]=a[3*j+h]-a[3*i+h];vv[h]=p[h]-a[3*i+h];}
    double len=dot3(ee2,ee2),s=len>0?dot3(vv,ee2)/len:0.;s=fmin(1.,fmax(0.,s));
    for(int h=0;h<3;h++){cp[h]=a[3*i+h]+s*ee2[h];d[h]=p[h]-cp[h];}
    double ds=dot3(d,d);
    if(ds<distance){distance=ds;for(int h=0;h<3;h++){best[h]=cp[h];bw[h]=0.;}bw[i]=1.-s;bw[j]=s;}
  }
}
__device__ double box_distance(const double* p,const double* box) {
  double d[3];for(int k=0;k<3;k++)d[k]=fmax(fmax(box[k]-p[k],p[k]-box[k+3]),0.);
  return sqrt(dot3(d,d));
}
__device__ double sphere_triangle(const double* p,double radius,const double* tri,const double* inward,
    double* wp,double* pp,double* n,double* bary) {
  triangle_closest(p,tri,wp,bary);double d[3];for(int k=0;k<3;k++)d[k]=p[k]-wp[k];
  double dist=sqrt(dot3(d,d)),gap=dist-radius;bool face=true;
  for(int k=0;k<3;k++){n[k]=dist>0?d[k]/dist:inward[k];if(bary[k]<=128.*DBL_EPSILON)face=false;}
  // Match capsule_triangle_many's witness before the final inward flip.
  if(face && ein3(d,inward)<0)for(int k=0;k<3;k++)n[k]=-inward[k];
  else if(face)for(int k=0;k<3;k++)n[k]=inward[k];
  for(int k=0;k<3;k++)pp[k]=face?p[k]-radius*n[k]:wp[k]+gap*n[k];
  if(ein3(n,inward)<0)for(int k=0;k<3;k++)n[k]=-n[k];
  if(face)for(int k=0;k<3;k++)n[k]=inward[k];
  return gap;
}
__global__ void wall_kernel(const double* x,const double* r,i64 n,const double* tri,const double* inward,
    const double* boxes,const i64* links,const i64* order,const i64* rim,double pad,
    i64* nearest,double* result,i64* contacts,double* contact_n,i64* counts) {
  i64 row=blockIdx.x*blockDim.x+threadIdx.x;if(row>=n)return;const double* p=x+3*row;
  double radius=r[row],best=INFINITY; i64 best_id=-1;int stack[128],top=0;stack[top++]=0;
  double bwp[3],bpp[3],bn[3],bb[3];
  while(top){int node=stack[--top];if(box_distance(p,boxes+6*node)>best+radius+pad)continue;
    const i64* l=links+4*node;
    if(l[0]>=0){double d0=box_distance(p,boxes+6*l[0]),d1=box_distance(p,boxes+6*l[1]);
      if(d0<d1){stack[top++]=l[1];stack[top++]=l[0];}else{stack[top++]=l[0];stack[top++]=l[1];}}
    else for(i64 z=l[2];z<l[2]+l[3];z++){
      i64 id=order[z];double wp[3],pp[3],normal[3],bary[3];double g=sphere_triangle(p,radius,tri+9*id,inward+3*id,wp,pp,normal,bary);
      if(g<best || (g==best && id<best_id)){best=g;best_id=id;
        for(int k=0;k<3;k++){bwp[k]=wp[k];bpp[k]=pp[k];bn[k]=normal[k];bb[k]=bary[k];}}
    }
  }
  nearest[row]=best_id;double* out=result+14*row;out[0]=best;
  int active=0,missing=-1,maximum=0;for(int k=0;k<3;k++){
    out[1+k]=bwp[k];out[4+k]=bpp[k];out[7+k]=bn[k];out[10+k]=bb[k];
    if(bb[k]>128.*DBL_EPSILON)active++;else missing=k;if(bb[k]>bb[maximum])maximum=k;
  }
  out[13]=active==2?rim[6*best_id+3+missing]:active<=1?rim[6*best_id+maximum]:0;
  double lower=fmax(2e-9,1e-3*radius);top=0;stack[top++]=0;int count=0;
  while(top){int node=stack[--top];if(box_distance(p,boxes+6*node)>radius+lower+pad)continue;
    const i64* l=links+4*node;if(l[0]>=0){stack[top++]=l[1];stack[top++]=l[0];}
    else for(i64 z=l[2];z<l[2]+l[3];z++){
      i64 id=order[z];double wp[3],pp[3],normal[3],bary[3];double g=sphere_triangle(p,radius,tri+9*id,inward+3*id,wp,pp,normal,bary);
      if(fabs(g-lower)<=pad){if(count<32){contacts[row*32+count]=id;for(int k=0;k<3;k++)contact_n[(row*32+count)*3+k]=normal[k];}count++;}
    }
  }
  // Canonical triangle order, independent of BVH construction/traversal order.
  int kept=min(count,32);for(int i=1;i<kept;i++){i64 id=contacts[row*32+i];double no[3];
    for(int k=0;k<3;k++)no[k]=contact_n[(row*32+i)*3+k];int j=i-1;
    while(j>=0&&contacts[row*32+j]>id){contacts[row*32+j+1]=contacts[row*32+j];
      for(int k=0;k<3;k++)contact_n[(row*32+j+1)*3+k]=contact_n[(row*32+j)*3+k];j--;}
    contacts[row*32+j+1]=id;for(int k=0;k<3;k++)contact_n[(row*32+j+1)*3+k]=no[k];}
  counts[row]=count;
}
std::vector<torch::Tensor> wall_query_cuda(std::vector<torch::Tensor> a,double roundoff) {
  TORCH_CHECK(a.size()==8,"Wall tensor bundle size");c10::cuda::CUDAGuard guard(a[0].device());
  i64 n=a[0].size(0);auto idopt=a[5].options();auto nearest=torch::empty({n},idopt);
  auto result=torch::empty({n,14},a[0].options());auto contacts=torch::full({n,32},-1,idopt);
  auto normals=torch::zeros({n,32,3},a[0].options());auto counts=torch::zeros({n},idopt);
  if(n)wall_kernel<<<(n+127)/128,128,0,at::cuda::getCurrentCUDAStream()>>>(
    a[0].data_ptr<double>(),a[1].data_ptr<double>(),n,a[2].data_ptr<double>(),a[3].data_ptr<double>(),
    a[4].data_ptr<double>(),a[5].data_ptr<i64>(),a[6].data_ptr<i64>(),a[7].data_ptr<i64>(),roundoff,
    nearest.data_ptr<i64>(),result.data_ptr<double>(),contacts.data_ptr<i64>(),normals.data_ptr<double>(),counts.data_ptr<i64>());
  C10_CUDA_KERNEL_LAUNCH_CHECK();return {nearest,result,contacts,normals,counts};
}

__device__ double path_gap(const double* p,const double* v,double t,const double* tri,double radius,double lower,double* derivative) {
  double center[3],wp[3],b[3],d[3];for(int k=0;k<3;k++)center[k]=p[k]+t*v[k];
  triangle_closest(center,tri,wp,b);for(int k=0;k<3;k++)d[k]=center[k]-wp[k];
  *derivative=dot3(d,v);return sqrt(dot3(d,d))-radius-lower;
}
__global__ void wall_path_kernel(const double* x,const double* velocity,const double* radii,const double* times,
    i64 count,const double* tri,const double* inward,const double* box,const i64* links,const i64* order,
    double pad,double* output) {
  i64 row=blockIdx.x*blockDim.x+threadIdx.x;if(row>=count)return;
  const double* p=x+3*row;const double* v=velocity+3*row;double radius=radii[row],dt=times[row];
  double lower=fmax(2e-9,1e-3*radius),delta[3],end[3],mid[3];
  for(int k=0;k<3;k++){delta[k]=dt*v[k];end[k]=p[k]+delta[k];mid[k]=.5*(p[k]+end[k]);}
  double distance=sqrt(dot3(delta,delta)),envelope=radius+lower+.5*distance;
  double first=dt,minbound=INFINITY;int root_id=-1,status=1,stack[128],top=0,visited=0;stack[top++]=0;
  while(top){int node=stack[--top];if(box_distance(mid,box+6*node)>envelope+pad)continue;
    const i64* l=links+4*node;if(l[0]>=0){stack[top++]=l[1];stack[top++]=l[0];continue;}
    for(i64 z=l[2];z<l[2]+l[3];z++){
      i64 id=order[z];const double* t=tri+9*id;visited++;
      double wp0[3],wp1[3],pp[3],n0[3],n1[3],b[3];
      double g0=sphere_triangle(p,radius,t,inward+3*id,wp0,pp,n0,b),g1=sphere_triangle(end,radius,t,inward+3*id,wp1,pp,n1,b);
      double d[3];for(int k=0;k<3;k++)d[k]=p[k]-wp0[k];double dist=sqrt(dot3(d,d));
      double direction=dist>0?ein3(d,delta)/dist:0.;
      // Same convex finite-triangle minimum and bisection root as CPU handoff.
      if(distance>0 && g0-lower+fmin(direction,0.) < -pad && g0-lower>pad){
        double deriv,eg=path_gap(p,v,dt,t,radius,lower,&deriv),upper=dt;
        if(eg>=-pad){
          if(deriv<=0)upper=0.;else{double lo=0.,hi=dt;
            for(int it=0;it<64;it++){double m=lo+(hi-lo)/2.;if(m==lo||m==hi)break;
              path_gap(p,v,m,t,radius,lower,&deriv);if(deriv<0)lo=m;else hi=m;}
            upper=lo+(hi-lo)/2.;}}
        double gupper=path_gap(p,v,upper,t,radius,lower,&deriv);
        if(gupper< -pad && g0-lower>pad){double lo=0.,hi=upper;
          for(int it=0;it<64;it++){double m=lo+(hi-lo)/2.;if(m==lo||m==hi)break;
            double g=path_gap(p,v,m,t,radius,lower,&deriv);if(g>=0)lo=m;else hi=m;
            if(g>=0 && g<=pad/4.){lo=m;break;}}
          double value=path_gap(p,v,lo,t,radius,lower,&deriv);
          if(value>=0 && value<=pad && (lo<first || (lo==first&&id<root_id))){first=lo;root_id=id;}}
      }
      if(fmin(g0,g1)-lower< -pad){status=-1;continue;}
      double bound=fmin(g0,g1)-distance-lower;
      if(bound < -pad){
        bound=-INFINITY;
        for(int nn=0;nn<2;nn++){const double* n=nn==0?n0:n1;double rel0[3],rel1[3],e1[3],e2[3];
          for(int k=0;k<3;k++){rel0[k]=p[k]-t[k];rel1[k]=end[k]-t[k];e1[k]=t[3+k]-t[k];e2[k]=t[6+k]-t[k];}
          double support=fmax(0.,fmax(dot3(e1,n),dot3(e2,n)));
          double h0=dot3(n,rel0)-radius-support,h1=dot3(n,rel1)-radius-support;
          bound=fmax(bound,fmin(h0,h1)-lower);}
        if(bound< -pad && status>0)status=0;
      }
      minbound=fmin(minbound,bound);
    }
  }
  double* out=output+5*row;out[0]=first;out[1]=root_id;out[2]=status;out[3]=minbound;out[4]=visited;
}
torch::Tensor wall_path_cuda(std::vector<torch::Tensor> a,double roundoff) {
  TORCH_CHECK(a.size()==10,"Path tensor bundle size");c10::cuda::CUDAGuard guard(a[0].device());i64 n=a[0].size(0);
  auto out=torch::empty({n,5},a[0].options());
  if(n)wall_path_kernel<<<(n+127)/128,128,0,at::cuda::getCurrentCUDAStream()>>>(a[0].data_ptr<double>(),
    a[1].data_ptr<double>(),a[2].data_ptr<double>(),a[3].data_ptr<double>(),n,a[4].data_ptr<double>(),a[5].data_ptr<double>(),
    a[6].data_ptr<double>(),a[7].data_ptr<i64>(),a[8].data_ptr<i64>(),roundoff,out.data_ptr<double>());
  C10_CUDA_KERNEL_LAUNCH_CHECK();return out;
}

#include "rigid_math.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numeric>
namespace rigid {
constexpr double mu=.001;const double pi=std::acos(-1.);
Vec add(Vec a,Vec b){for(int k=0;k<3;k++)a[k]+=b[k];return a;}Vec sub(Vec a,Vec b){return add(a,mul(b,-1));}Vec mul(Vec a,double s){for(double&v:a)v*=s;return a;}double dot(Vec a,Vec b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}double norm(Vec a){return std::sqrt(dot(a,a));}Vec cross(Vec a,Vec b){return{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}Vec mv(Mat a,Vec x){Vec r{};for(int k=0;k<3;k++)for(int l=0;l<3;l++)r[k]+=a[3*k+l]*x[l];return r;}
Background FlowGradientSampler::query(Vec x)const{
 auto old=frozen::FlowFieldSampler(field).query(x);if(old.status!=frozen::Status::VALID)throw std::runtime_error("INVALID_FLOW_QUERY");Background b;b.u=old.u;
 Vec t;long long base[3];for(int k=0;k<3;k++){double q=(x[k]-field.origin[k])/field.dx;base[k]=std::floor(q);t[k]=q-base[k];}
 // Analytic derivative, subtract one common corner velocity to suppress cancellation.
 Vec anchor{};size_t anchor_slot;auto st=field.node(base[0]+field.dims[0]*(base[1]+field.dims[1]*base[2]),anchor_slot);if(st!=frozen::Status::VALID)throw std::runtime_error("INVALID_GRADIENT_CORNER");for(int k=0;k<3;k++)anchor[k]=field.velocity[3*anchor_slot+k];
 for(int z=0;z<2;z++)for(int y=0;y<2;y++)for(int xb=0;xb<2;xb++){
  int bit[3]={xb,y,z};uint64_t id=base[0]+xb+field.dims[0]*(base[1]+y+field.dims[1]*(base[2]+z));size_t slot;st=field.node(id,slot);if(st!=frozen::Status::VALID)throw std::runtime_error("INVALID_GRADIENT_CORNER");
  for(int axis=0;axis<3;axis++){double weight=(bit[axis]?1.:-1.)/field.dx;for(int k=0;k<3;k++)if(k!=axis)weight*=bit[k]?t[k]:1-t[k];for(int k=0;k<3;k++)b.grad[3*k+axis]+=weight*(field.velocity[3*slot+k]-anchor[k]);}
 }
 for(int k=0;k<3;k++)for(int l=0;l<3;l++)b.strain[3*k+l]=.5*(b.grad[3*k+l]+b.grad[3*l+k]);b.omega={.5*(b.grad[7]-b.grad[5]),.5*(b.grad[2]-b.grad[6]),.5*(b.grad[3]-b.grad[1])};return b;
}
Coeff coefficients(double a,double b,double h){
 if(!(a>0&&b>0&&h>=-1e-12))throw std::runtime_error("INVALID_PAIR_GEOMETRY");
 const double s=a+b,re=a*b/s,he=std::max(h,.001*re),L=std::log(s/(2*he)),ab=a*b;
 // Symmetric dimensional polynomials, independently coded from beta/xi Python reference.
 double sq=6*pi*mu*(ab*ab/(s*s*he)+ab*(a*a+7*ab+b*b)/(5*s*s*s)*L+(std::pow(a,4)+18*a*a*ab-29*ab*ab+18*b*b*ab+std::pow(b,4))/(21*std::pow(s,4))*he*L);
 double sh0=8*pi*mu*ab*(2*a*a+ab+2*b*b)/(5*s*s*s)*L;
 double sh=sh0+8*pi*mu*(16*std::pow(a,4)-45*a*a*ab+58*ab*ab-45*b*b*ab+16*std::pow(b,4))/(125*std::pow(s,4))*he*L;
 // Effective torque arms from the independently published YB/YA ratio.
 double bi=a*(4*a+b),bj=b*(a+4*b),li=(s+h)*bi/(bi+bj),lj=(s+h)*bj/(bi+bj);
 // Schur complement C11 - B11^2 / YA, simplified symbolically.
 double pu=6*pi*mu*a*a*a*b*b*b/(s*(2*a*a+ab+2*b*b))*L*std::pow((s+h)/s,2);
 return{sq,sh,pu,0,li,lj};
}
static void rankadd(Array&R,Array&rhs,const Array&v,double w,double target){int n=v.size();for(int k=0;k<n;k++){rhs[k]+=w*v[k]*target;for(int l=0;l<n;l++)R[k*n+l]+=w*v[k]*v[l];}}
System assemble(const std::vector<Particle>&p,const std::vector<Background>&bg,const std::vector<std::pair<int,int>>&pairs){
 int d=6*p.size();System s;s.R.assign(d*d,0);s.rhs.assign(d,0);Array q0(d);
 for(size_t i=0;i<p.size();i++)for(int k=0;k<6;k++){double drag=k<3?6*pi*mu*p[i].a:8*pi*mu*std::pow(p[i].a,3);int l=6*i+k;q0[l]=k<3?bg[i].u[k]:bg[i].omega[k-3];s.R[l*d+l]=drag;s.rhs[l]=drag*(q0[l]+(k<3?bg[i].du[k]:bg[i].dw[k-3]));}
 for(auto ij:pairs){int i=ij.first,j=ij.second;Vec rvec=sub(p[i].x,p[j].x);double r=norm(rvec),h=r-p[i].a-p[j].a;double cutoff=.2*p[i].a*p[j].a/(p[i].a+p[j].a);if(h>=cutoff)continue;Vec n=mul(rvec,1/r);Coeff c=coefficients(p[i].a,p[j].a,h);Mat E;for(int k=0;k<9;k++)E[k]=.5*(bg[i].strain[k]+bg[j].strain[k]);Vec strain=mul(mv(E,n),p[i].a+p[j].a);
  for(int k=0;k<3;k++){
   Array b(d),pump(d);for(int l=0;l<3;l++){double P=(k==l?1.:0.)-n[k]*n[l];b[6*i+l]=P;b[6*j+l]=-P;Vec unit{};unit[l]=1;double crossnl=cross(n,unit)[k];b[6*i+3+l]=c.li*crossnl;b[6*j+3+l]=c.lj*crossnl;pump[6*i+3+l]=P;pump[6*j+3+l]=-P;}
   double target=std::inner_product(b.begin(),b.end(),q0.begin(),0.)-(strain[k]-n[k]*dot(n,strain));rankadd(s.R,s.rhs,b,c.sh,target);rankadd(s.R,s.rhs,pump,c.pu,std::inner_product(pump.begin(),pump.end(),q0.begin(),0.));
  }
  Array b(d);for(int k=0;k<3;k++){b[6*i+k]=n[k];b[6*j+k]=-n[k];}rankadd(s.R,s.rhs,b,c.sq,std::inner_product(b.begin(),b.end(),q0.begin(),0.)-dot(n,strain));Event e{};e.i=i;e.j=j;e.gap=h;e.n=n;e.c=c;s.events.push_back(e);
 }
 return s;
}
static double dp(const Array&a,const Array&b){return std::inner_product(a.begin(),a.end(),b.begin(),0.);}static Array matvec(const Array&A,const Array&x){int n=x.size();Array r(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)r[i]+=A[i*n+j]*x[j];return r;}
static Array pcg(const Array&A,const Array&b,int&iters){int n=b.size();Array x(n),r=b,z(n),v;double bn=std::sqrt(dp(b,b));if(bn==0)return x;for(int i=0;i<n;i++)z[i]=r[i]/A[i*n+i];v=z;double rz=dp(r,z);
 for(int it=0;it<10*n+100;it++){auto av=matvec(A,v);double curvature=dp(v,av);if(!(curvature>0))throw std::runtime_error("NON_POSITIVE_PCG_CURVATURE");double alpha=rz/curvature;for(int i=0;i<n;i++){x[i]+=alpha*v[i];r[i]-=alpha*av[i];}iters++;if(std::sqrt(dp(r,r))/bn<2e-14)return x;for(int i=0;i<n;i++)z[i]=r[i]/A[i*n+i];double next=dp(r,z),beta=next/rz;for(int i=0;i<n;i++)v[i]=z[i]+beta*v[i];rz=next;}
 throw std::runtime_error("PCG_DID_NOT_CONVERGE");}
static Array dense(Array A,Array b){int n=b.size();for(int k=0;k<n;k++){int pivot=k;for(int j=k+1;j<n;j++)if(std::abs(A[j*n+k])>std::abs(A[pivot*n+k]))pivot=j;if(std::abs(A[pivot*n+k])<1e-30)throw std::runtime_error("SINGULAR_CONSTRAINT_SYSTEM");for(int j=k;j<n;j++)std::swap(A[k*n+j],A[pivot*n+j]);std::swap(b[k],b[pivot]);double a=A[k*n+k];for(int j=k;j<n;j++)A[k*n+j]/=a;b[k]/=a;for(int i=0;i<n;i++)if(i!=k){double t=A[i*n+k];for(int j=k;j<n;j++)A[i*n+j]-=t*A[k*n+j];b[i]-=t*b[k];}}return b;}
Solution solve(const std::vector<Particle>&p,const std::vector<Background>&bg,const std::vector<std::pair<int,int>>&pairs,const std::vector<Particle>*base,double dt,const std::vector<WallBlock>*walls){
 auto sys=assemble(p,bg,pairs);int d=6*p.size();
 // Add ONLY wall excess and its resistance-only ambient-slip right-hand side.
 if(walls)for(auto const&w:*walls){int i=w.particle_index;for(int k=0;k<6;k++)for(int l=0;l<6;l++){
  double r=w.excess[6*k+l];sys.R[(6*i+k)*d+6*i+l]+=r;
  sys.rhs[6*i+k]+=r*(l<3?bg[i].u[l]:bg[i].omega[l-3]);}}
Array scale(d),A(d*d),b(d);for(int k=0;k<d;k++)scale[k]=k%6<3?1:1/p[k/6].a;for(int k=0;k<d;k++){b[k]=sys.rhs[k]*scale[k]/1e-8;for(int l=0;l<d;l++)A[k*d+l]=sys.R[k*d+l]*scale[k]*scale[l]/1e-8;}Solution out;Array y=pcg(A,b,out.iterations),original=y;std::vector<Array> J,Z;Array bounds,lambda;std::vector<int>pairids;
 if(base&&dt>0){
  // Unified active set. Wall rows first in (particle ID, triangle ID) order.
  // Pair rows retain the frozen (minID,maxID) order. Responses are computed
  // lazily only when a row activates; this does not change the equations.
  if(walls)for(auto const&w:*walls)if(w.hard){Array row(d);for(int l=0;l<3;l++)row[6*w.particle_index+l]=w.base_normal[l];J.push_back(row);Z.push_back({});bounds.push_back(-std::max(w.base_gap,0.)/dt);pairids.push_back(-1-w.particle_index);}
  for(size_t k=0;k<pairs.size();k++){int i=pairs[k].first,j=pairs[k].second;Vec rv=sub((*base)[i].x,(*base)[j].x);double r=norm(rv),h=r-p[i].a-p[j].a;if(h>=.2*p[i].a*p[j].a/(p[i].a+p[j].a))continue;Vec n=mul(rv,1/r);Array row(d);for(int l=0;l<3;l++){row[6*i+l]=n[l];row[6*j+l]=-n[l];}J.push_back(row);Z.push_back({});bounds.push_back(-std::max(h,0.)/dt);pairids.push_back(k);}
  for(int iter=0;iter<4*int(J.size()+1)+20;iter++){
   int worst=-1;double v=-1e-14;for(size_t k=0;k<J.size();k++){double e=dp(J[k],y)-bounds[k];if(e<v&&std::find(out.constraints.begin(),out.constraints.end(),k)==out.constraints.end()){worst=k;v=e;}}if(worst<0)break;out.constraints.push_back(worst);if(Z[worst].empty())Z[worst]=pcg(A,J[worst],out.iterations);
   for(int inner=0;inner<100;inner++){int m=out.constraints.size();Array G(m*m),rhs(m);for(int k=0;k<m;k++){int ik=out.constraints[k];rhs[k]=bounds[ik]-dp(J[ik],original);for(int l=0;l<m;l++)G[k*m+l]=dp(J[ik],Z[out.constraints[l]]);}lambda=dense(G,rhs);int neg=-1;double value=-1e-14;for(int k=0;k<m;k++)if(lambda[k]<value){value=lambda[k];neg=k;}if(neg>=0){out.constraints.erase(out.constraints.begin()+neg);continue;}y=original;for(int k=0;k<m;k++)for(int l=0;l<d;l++)y[l]+=lambda[k]*Z[out.constraints[k]][l];break;}
  }
  for(size_t k=0;k<J.size();k++)if(dp(J[k],y)<bounds[k]-1e-12)throw std::runtime_error("CONSTRAINT_VIOLATION");
 }
 auto res=matvec(A,y);for(int k=0;k<d;k++)res[k]-=b[k];for(size_t k=0;k<out.constraints.size();k++)for(int l=0;l<d;l++)res[l]-=lambda[k]*J[out.constraints[k]][l];out.residual=std::sqrt(dp(res,res))/std::max(std::sqrt(dp(b,b)),1e-30);if(out.residual>1e-10)throw std::runtime_error("RESIDUAL_GATE_FAILED");out.q.resize(d);for(int k=0;k<d;k++)out.q[k]=y[k]*scale[k];out.events=sys.events;
 for(auto&e:out.events){Vec v{},w{},wi{},wj{};for(int k=0;k<3;k++){v[k]=out.q[6*e.i+k]-out.q[6*e.j+k]-bg[e.i].u[k]+bg[e.j].u[k];wi[k]=out.q[6*e.i+3+k]-bg[e.i].omega[k];wj[k]=out.q[6*e.j+3+k]-bg[e.j].omega[k];}Mat E;for(int k=0;k<9;k++)E[k]=.5*(bg[e.i].strain[k]+bg[e.j].strain[k]);v=add(v,mul(mv(E,e.n),p[e.i].a+p[e.j].a));w=sub(wi,wj);e.normal=mul(e.n,dot(e.n,v));e.slip=add(sub(v,e.normal),cross(e.n,add(mul(wi,e.c.li),mul(wj,e.c.lj))));e.transverse=sub(w,mul(e.n,dot(e.n,w)));e.fn=mul(e.normal,-e.c.sq);e.fs=mul(e.slip,-e.c.sh);e.force=add(e.fn,e.fs);e.tpu=mul(e.transverse,-e.c.pu);e.ti=add(mul(cross(e.n,e.fs),-e.c.li),e.tpu);e.tj=sub(mul(cross(e.n,e.fs),-e.c.lj),e.tpu);for(int k:out.constraints){if(pairids[k]<0)continue;auto ij=pairs[pairids[k]];if(ij.first==e.i&&ij.second==e.j)e.constraint=true;}}
 for(int k:out.constraints)if(pairids[k]<0)out.wall_constraints.push_back(-1-pairids[k]);
 return out;
}
}
extern "C" void rigid_coefficients(int n,const double*in,double*out){for(int i=0;i<n;i++){auto c=rigid::coefficients(in[3*i],in[3*i+1],in[3*i+2]);out[6*i]=c.sq;out[6*i+1]=c.sh;out[6*i+2]=c.pu;out[6*i+3]=c.tw;out[6*i+4]=c.li;out[6*i+5]=c.lj;}}
extern "C" int rigid_audit_system(int n,const double*state,int np,const int*pair,double*matrix,double*rhs,double*solution,double*stats){try{std::vector<rigid::Particle>p(n);std::vector<rigid::Background>b(n);for(int i=0;i<n;i++){p[i].id=i;p[i].a=state[19*i+3];for(int k=0;k<3;k++){p[i].x[k]=state[19*i+k];b[i].u[k]=state[19*i+4+k];b[i].omega[k]=state[19*i+7+k];}for(int k=0;k<9;k++)b[i].strain[k]=state[19*i+10+k];}std::vector<std::pair<int,int>>pairs;for(int k=0;k<np;k++)pairs.emplace_back(pair[2*k],pair[2*k+1]);auto s=rigid::assemble(p,b,pairs);std::copy(s.R.begin(),s.R.end(),matrix);std::copy(s.rhs.begin(),s.rhs.end(),rhs);if(solution){auto q=rigid::solve(p,b,pairs);std::copy(q.q.begin(),q.q.end(),solution);stats[0]=q.residual;stats[1]=q.iterations;}return 0;}catch(...){return 1;}}
extern "C" int gradient_audit(const char*path,int n,const double*x,double*out){try{frozen::FlowField f(path);rigid::FlowGradientSampler s(f);for(int i=0;i<n;i++){auto b=s.query({x[3*i],x[3*i+1],x[3*i+2]});for(int k=0;k<3;k++){out[15*i+k]=b.u[k];out[15*i+3+k]=b.omega[k];}for(int k=0;k<9;k++)out[15*i+6+k]=b.grad[k];}return 0;}catch(...){return 1;}}

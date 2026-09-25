#include "wall_distance.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <limits>
namespace passive {
Vec minus(Vec a,Vec b){for(int j=0;j<3;++j)a[j]-=b[j];return a;}
Vec plus_scaled(Vec a,Vec b,double s){for(int j=0;j<3;++j)a[j]+=s*b[j];return a;}
static double dot(Vec a,Vec b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm(Vec a){return std::sqrt(dot(a,a));}
static double edge2(Vec p,Vec a,Vec b){Vec v=minus(b,a),w=minus(p,a);double vv=dot(v,v);double t=vv>0?std::clamp(dot(w,v)/vv,0.,1.):0.;Vec r=minus(w,plus_scaled(Vec{0,0,0},v,t));return dot(r,r);}
static double triangle2(Vec p,Triangle const&t){
 Vec ab=minus(t.b,t.a),ac=minus(t.c,t.a),ap=minus(p,t.a);double aa=dot(ab,ab),bb=dot(ab,ac),cc=dot(ac,ac),pa=dot(ap,ab),pc=dot(ap,ac),det=aa*cc-bb*bb;
 if(det>0){double v=(cc*pa-bb*pc)/det,w=(aa*pc-bb*pa)/det;if(v>=0&&w>=0&&v+w<=1){Vec q=plus_scaled(plus_scaled(t.a,ab,v),ac,w);Vec e=minus(p,q);return dot(e,e);}}
 return std::min({edge2(p,t.a,t.b),edge2(p,t.b,t.c),edge2(p,t.c,t.a)});
}
WallDistance::WallDistance(const std::string& path){std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("Cannot open closed binary STL");char header[80];uint32_t n=0;in.read(header,80);in.read(reinterpret_cast<char*>(&n),4);if(!in||!n||n>10000000)throw std::runtime_error("Invalid binary STL");triangles.reserve(n);
 for(uint32_t i=0;i<n;++i){float row[12];uint16_t attribute;in.read(reinterpret_cast<char*>(row),48);in.read(reinterpret_cast<char*>(&attribute),2);if(!in)throw std::runtime_error("Truncated STL");Triangle t;for(int j=0;j<3;++j){t.a[j]=row[3+j];t.b[j]=row[6+j];t.c[j]=row[9+j];t.center[j]=(t.a[j]+t.b[j]+t.c[j])/3;}triangles.push_back(t);}build(0,n);}
int WallDistance::build(int b,int e){Node n;n.begin=b;n.end=e;n.lo={INFINITY,INFINITY,INFINITY};n.hi={-INFINITY,-INFINITY,-INFINITY};for(int i=b;i<e;++i)for(auto const&v:{triangles[i].a,triangles[i].b,triangles[i].c})for(int j=0;j<3;++j){n.lo[j]=std::min(n.lo[j],v[j]);n.hi[j]=std::max(n.hi[j],v[j]);}int idx=nodes.size();nodes.push_back(n);if(e-b>8){int axis=0;for(int j=1;j<3;++j)if(n.hi[j]-n.lo[j]>n.hi[axis]-n.lo[axis])axis=j;int mid=(b+e)/2;std::nth_element(triangles.begin()+b,triangles.begin()+mid,triangles.begin()+e,[axis](auto const&a,auto const&b){return a.center[axis]<b.center[axis];});int l=build(b,mid),r=build(mid,e);nodes[idx].left=l;nodes[idx].right=r;}return idx;}
static double box2(Node const&n,Vec p){double d=0;for(int j=0;j<3;++j){double v=p[j]<n.lo[j]?n.lo[j]-p[j]:(p[j]>n.hi[j]?p[j]-n.hi[j]:0);d+=v*v;}return d;}
void WallDistance::search(int i,Vec p,double&best) const{auto const&n=nodes[i];if(box2(n,p)>best)return;if(n.left<0){for(int k=n.begin;k<n.end;++k)best=std::min(best,triangle2(p,triangles[k]));return;}int a=n.left,b=n.right;if(box2(nodes[b],p)<box2(nodes[a],p))std::swap(a,b);search(a,p,best);search(b,p,best);}
double WallDistance::distance(Vec p)const{double best=INFINITY;search(0,p,best);return std::sqrt(best);}
bool WallDistance::segment_safe(Vec a,Vec b,double minimum,int depth)const{double da=distance(a),db=distance(b),len=norm(minus(b,a));if(std::min(da,db)<minimum)return false;if(std::min(da,db)-len/2>=minimum)return true;if(depth>=20)return false;Vec mid=plus_scaled(a,minus(b,a),.5);return segment_safe(a,mid,minimum,depth+1)&&segment_safe(mid,b,minimum,depth+1);}
}

namespace passive {
static Vec closest_edge(Vec p,Vec a,Vec b){Vec u=minus(b,a);double den=dot(u,u);return plus_scaled(a,u,den>0?std::clamp(dot(minus(p,a),u)/den,0.,1.):0.);}
static Vec closest_triangle(Vec p,Triangle const&t){
 Vec ab=minus(t.b,t.a),ac=minus(t.c,t.a),ap=minus(p,t.a);double aa=dot(ab,ab),bb=dot(ab,ac),cc=dot(ac,ac),pa=dot(ap,ab),pc=dot(ap,ac),det=aa*cc-bb*bb;
 if(det>0){double v=(cc*pa-bb*pc)/det,w=(aa*pc-bb*pa)/det;if(v>=0&&w>=0&&v+w<=1)return plus_scaled(plus_scaled(t.a,ab,v),ac,w);}
 Vec best=closest_edge(p,t.a,t.b);for(auto q:{closest_edge(p,t.b,t.c),closest_edge(p,t.c,t.a)})if(norm(minus(p,q))<norm(minus(p,best)))best=q;return best;
}
void WallDistance::closest_search(int i,Vec p,ClosestWallPoint&best) const {
 auto const&n=nodes[i];double tol=64*std::numeric_limits<double>::epsilon()*std::max(norm(p),1e-6);
 if(box2(n,p)>(best.distance+tol)*(best.distance+tol))return;
 if(n.left<0){for(int k=n.begin;k<n.end;k++){Vec q=closest_triangle(p,triangles[k]);double d=norm(minus(p,q));
  if(d<best.distance-tol){best={q,d,false};}
  else if(std::abs(d-best.distance)<=tol){if(norm(minus(q,best.point))>8*tol)best.ambiguous=true;if(d<best.distance){best.point=q;best.distance=d;}}
 }return;}
 int a=n.left,b=n.right;if(box2(nodes[b],p)<box2(nodes[a],p))std::swap(a,b);closest_search(a,p,best);closest_search(b,p,best);
}
ClosestWallPoint WallDistance::closest(Vec p)const{ClosestWallPoint best;best.distance=INFINITY;closest_search(0,p,best);return best;}
static double segment_segment2(Vec p,Vec q,Vec a,Vec b){
 Vec u=minus(q,p),v=minus(b,a),w=minus(p,a);double A=dot(u,u),B=dot(u,v),C=dot(v,v),D=dot(u,w),E=dot(v,w),det=A*C-B*B;
 double best=std::min({edge2(p,a,b),edge2(q,a,b),edge2(a,p,q),edge2(b,p,q)});
 if(det>0){double s=(B*E-C*D)/det,t=(A*E-B*D)/det;if(s>=0&&s<=1&&t>=0&&t<=1){Vec d=minus(plus_scaled(p,u,s),plus_scaled(a,v,t));best=std::min(best,dot(d,d));}}
 return best;
}
static Vec cross(Vec a,Vec b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
static double segment_triangle2(Vec a,Vec b,Triangle const&t){
 Vec d=minus(b,a),n=cross(minus(t.b,t.a),minus(t.c,t.a));double den=dot(d,n);
 if(den!=0){double f=dot(minus(t.a,a),n)/den;if(f>=0&&f<=1){Vec x=plus_scaled(a,d,f),v=minus(t.b,t.a),w=minus(t.c,t.a),h=minus(x,t.a);double vv=dot(v,v),vw=dot(v,w),ww=dot(w,w),hv=dot(h,v),hw=dot(h,w),det=vv*ww-vw*vw;if(det>0){double u=(ww*hv-vw*hw)/det,z=(vv*hw-vw*hv)/det;if(u>=0&&z>=0&&u+z<=1)return 0;}}}
 return std::min({triangle2(a,t),triangle2(b,t),segment_segment2(a,b,t.a,t.b),segment_segment2(a,b,t.b,t.c),segment_segment2(a,b,t.c,t.a)});
}
bool WallDistance::segment_search(int i,Vec a,Vec b,double minimum)const{
 auto const&n=nodes[i];double lower2=0;for(int j=0;j<3;j++){double lo=std::min(a[j],b[j]),hi=std::max(a[j],b[j]),gap=std::max({n.lo[j]-hi,lo-n.hi[j],0.});lower2+=gap*gap;}if(lower2>=minimum*minimum)return true;
 if(n.left<0){for(int k=n.begin;k<n.end;k++)if(segment_triangle2(a,b,triangles[k])<minimum*minimum)return false;return true;}
 return segment_search(n.left,a,b,minimum)&&segment_search(n.right,a,b,minimum);
}
bool WallDistance::segment_safe_exact(Vec a,Vec b,double minimum)const{return segment_search(0,a,b,minimum);}
}

#include "wall_distance.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <Eigen/Dense>
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
 for(uint32_t i=0;i<n;++i){float row[12];uint16_t attribute;in.read(reinterpret_cast<char*>(row),48);in.read(reinterpret_cast<char*>(&attribute),2);if(!in)throw std::runtime_error("Truncated STL");Triangle t;for(int j=0;j<3;++j){t.a[j]=row[3+j];t.b[j]=row[6+j];t.c[j]=row[9+j];t.center[j]=(t.a[j]+t.b[j]+t.c[j])/3;}t.original_id=i;Vec u=minus(t.b,t.a),v=minus(t.c,t.a);t.normal={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};double nn=norm(t.normal);if(nn<=0)throw std::runtime_error("DEGENERATE_STL_TRIANGLE");for(double &x:t.normal)x/=nn;triangles.push_back(t);}build(0,n);}
int WallDistance::build(int b,int e){Node n;n.begin=b;n.end=e;n.lo={INFINITY,INFINITY,INFINITY};n.hi={-INFINITY,-INFINITY,-INFINITY};for(int i=b;i<e;++i)for(auto const&v:{triangles[i].a,triangles[i].b,triangles[i].c})for(int j=0;j<3;++j){n.lo[j]=std::min(n.lo[j],v[j]);n.hi[j]=std::max(n.hi[j],v[j]);}int idx=nodes.size();nodes.push_back(n);if(e-b>8){int axis=0;for(int j=1;j<3;++j)if(n.hi[j]-n.lo[j]>n.hi[axis]-n.lo[axis])axis=j;int mid=(b+e)/2;std::nth_element(triangles.begin()+b,triangles.begin()+mid,triangles.begin()+e,[axis](auto const&a,auto const&b){return a.center[axis]!=b.center[axis]?a.center[axis]<b.center[axis]:a.original_id<b.original_id;});int l=build(b,mid),r=build(mid,e);nodes[idx].left=l;nodes[idx].right=r;}return idx;}
static double box2(Node const&n,Vec p){double d=0;for(int j=0;j<3;++j){double v=p[j]<n.lo[j]?n.lo[j]-p[j]:(p[j]>n.hi[j]?p[j]-n.hi[j]:0);d+=v*v;}return d;}
void WallDistance::search(int i,Vec p,double&best) const{auto const&n=nodes[i];if(box2(n,p)>best)return;if(n.left<0){for(int k=n.begin;k<n.end;++k)best=std::min(best,triangle2(p,triangles[k]));return;}int a=n.left,b=n.right;if(box2(nodes[b],p)<box2(nodes[a],p))std::swap(a,b);search(a,p,best);search(b,p,best);}
double WallDistance::distance(Vec p)const{double best=INFINITY;search(0,p,best);return std::sqrt(best);}
bool WallDistance::segment_safe(Vec a,Vec b,double minimum,int depth)const{double da=distance(a),db=distance(b),len=norm(minus(b,a));if(std::min(da,db)<minimum)return false;if(std::min(da,db)-len/2>=minimum)return true;if(depth>=30)return false;Vec mid=plus_scaled(a,minus(b,a),.5);return segment_safe(a,mid,minimum,depth+1)&&segment_safe(mid,b,minimum,depth+1);}
static Vec closest_edge(Vec p,Vec a,Vec b){Vec v=minus(b,a);double t=std::clamp(dot(minus(p,a),v)/dot(v,v),0.,1.);return plus_scaled(a,v,t);}
static Vec closest_triangle(Vec p,Triangle const&t){
 Vec ab=minus(t.b,t.a),ac=minus(t.c,t.a),ap=minus(p,t.a);double aa=dot(ab,ab),bb=dot(ab,ac),cc=dot(ac,ac),pa=dot(ap,ab),pc=dot(ap,ac),det=aa*cc-bb*bb;
 if(det>0){double v=(cc*pa-bb*pc)/det,w=(aa*pc-bb*pa)/det;if(v>=0&&w>=0&&v+w<=1)return plus_scaled(plus_scaled(t.a,ab,v),ac,w);}
 Vec best=closest_edge(p,t.a,t.b);double best2=dot(minus(p,best),minus(p,best));
 for(auto edge:{std::pair<Vec,Vec>{t.b,t.c},{t.c,t.a}}){Vec q=closest_edge(p,edge.first,edge.second);double d2=dot(minus(p,q),minus(p,q));if(d2<best2){best2=d2;best=q;}}
 return best;
}
void WallDistance::nearest(int i,Vec p,double&best,int&which,Vec&point)const{
 const auto &n=nodes[i];if(box2(n,p)>best)return;
 if(n.left<0){for(int k=n.begin;k<n.end;k++){Vec q=closest_triangle(p,triangles[k]);double d=dot(minus(p,q),minus(p,q));if(d<best||(d==best&&(which<0||triangles[k].original_id<triangles[which].original_id))){best=d;which=k;point=q;}}return;}
 int a=n.left,b=n.right;if(box2(nodes[b],p)<box2(nodes[a],p))std::swap(a,b);nearest(a,p,best,which,point);nearest(b,p,best,which,point);
}
void WallDistance::patch(int i,Vec p,double r2,std::vector<int>&out)const{
 const auto &n=nodes[i];if(box2(n,p)>r2)return;
 if(n.left<0){for(int k=n.begin;k<n.end;k++)if(triangle2(p,triangles[k])<=r2)out.push_back(k);return;}
 patch(n.left,p,r2,out);patch(n.right,p,r2,out);
}
WallQuery WallDistance::query(Vec p,double radius)const{
 WallQuery q;if(!(radius>0))return q;for(double x:p)if(!std::isfinite(x))return q;
 double best=INFINITY;int k=-1;nearest(0,p,best,k,q.closest);if(k<0||!(best>0))return q;
 q.distance=std::sqrt(best);q.gap=q.distance-radius;q.triangle_id=triangles[k].original_id;q.triangle_normal=triangles[k].normal;
 q.normal=minus(p,q.closest);for(double&v:q.normal)v/=q.distance;
 double projection=dot(q.normal,q.triangle_normal);q.triangle_normal_fluid=q.triangle_normal; if(projection<0)for(double&v:q.triangle_normal_fluid)v=-v;
 constexpr double deg=180./3.14159265358979323846;q.nearest_normal_angle_deg=std::acos(std::clamp(std::abs(projection),0.,1.))*deg;q.ambiguous=q.nearest_normal_angle_deg>20.;
 std::vector<int>patch_ids;patch(0,q.closest,4*radius*radius,patch_ids);
 std::sort(patch_ids.begin(),patch_ids.end(),[&](int i,int j){return triangles[i].original_id<triangles[j].original_id;});q.patch_triangles=patch_ids.size();
 // Deterministic area-weighted PCA of whole triangles intersecting the 2a ball.
 // Each triangle supplies three vertex samples with weight area/3. This exact
 // discretized patch rule is repeated independently in the Python finalizer.
 Eigen::Vector3d mean=Eigen::Vector3d::Zero();Eigen::Matrix3d second=Eigen::Matrix3d::Zero();double weight=0;
 std::vector<std::pair<Vec,double>>normals;
 for(int id:patch_ids){auto const&t=triangles[id];Vec u=minus(t.b,t.a),v=minus(t.c,t.a);Vec c={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};double area=.5*norm(c);normals.push_back({t.normal,area});
  for(Vec x:{t.a,t.b,t.c}){Eigen::Vector3d z;for(int j=0;j<3;j++)z[j]=(x[j]-q.closest[j])/radius;double w=area/3;mean+=w*z;second+=w*z*z.transpose();weight+=w;}}
 if(!(weight>0))return q;mean/=weight;Eigen::Matrix3d cov=second/weight-mean*mean.transpose();Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>eig(cov);
 if(eig.info()!=Eigen::Success||eig.eigenvalues()[1]<=1e-20)return q;
 Eigen::Vector3d pn=eig.eigenvectors().col(0);if(pn.dot(Eigen::Vector3d(q.normal[0],q.normal[1],q.normal[2]))<0)pn=-pn;
 for(int j=0;j<3;j++)q.plane_normal[j]=pn[j];q.rms_over_a=std::sqrt(std::max(0.,eig.eigenvalues()[0]));
 std::vector<std::pair<double,double>>angles;double total=0;for(auto const&n:normals){double a=std::acos(std::clamp(std::abs(dot(n.first,q.plane_normal)),0.,1.))*deg;angles.push_back({a,n.second});total+=n.second;}
 std::sort(angles.begin(),angles.end());double sum=0;for(auto aw:angles){sum+=aw.second;if(sum>=.95*total){q.normal_spread_deg=aw.first;break;}}
 q.valid=true;q.planar=q.rms_over_a<=.10&&q.normal_spread_deg<=20.&&!q.ambiguous;return q;
}
}

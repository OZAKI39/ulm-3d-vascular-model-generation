#pragma once
#include "rigid_math.hpp"
#include "wall_distance.hpp"
#include "adaptive_injection.hpp"
#include <fstream>
#include <memory>
namespace workflow {
using rigid::Vec;using rigid::sub;using rigid::add;using rigid::mul;using rigid::dot;using rigid::norm;using rigid::cross;
struct Tri {Vec a,b,c;};
inline std::vector<Tri> stl_triangles(std::string path){
 std::ifstream f(path,std::ios::binary);char head[80];uint32_t n;f.read(head,80);f.read(reinterpret_cast<char*>(&n),4);if(!f||n>10000000)throw std::runtime_error("INVALID_PORT_STL");std::vector<Tri>out;
 for(uint32_t i=0;i<n;i++){float x[12];uint16_t att;f.read(reinterpret_cast<char*>(x),48);f.read(reinterpret_cast<char*>(&att),2);Tri t;for(int k=0;k<3;k++){t.a[k]=x[3+k];t.b[k]=x[6+k];t.c[k]=x[9+k];}out.push_back(t);}if(!f)throw std::runtime_error("TRUNCATED_PORT_STL");return out;
}
struct Port {std::string name;Vec center,normal;std::vector<Tri>tri;};
struct Crossing {int port=-1;double fraction=2;Vec x{};};
class Geometry {
 public:passive::WallDistance wall;std::vector<std::unique_ptr<passive::WallDistance>>regions;std::vector<Port>ports;
 explicit Geometry(std::string root):wall(root+"/geometry/WALL_ONLY.stl"){
  for(int k=0;k<3;k++)regions.push_back(std::make_unique<passive::WallDistance>(root+"/geometry/REGION_"+std::to_string(k)+".stl"));
  std::ifstream f(root+"/geometry/ports.txt");Port p;std::string path;
  while(f>>p.name>>p.center[0]>>p.center[1]>>p.center[2]>>p.normal[0]>>p.normal[1]>>p.normal[2]>>path){p.tri=stl_triangles(root+"/"+path);ports.push_back(p);}
  if(ports.size()!=4)throw std::runtime_error("PORT_MANIFEST_COUNT");
 }
 int region(Vec x)const{int best=0;double d=regions[0]->distance(x);for(int k=1;k<3;k++){double v=regions[k]->distance(x);if(v<d){best=k;d=v;}}return best;}
 Crossing crossing(Vec start,Vec end)const{
  Crossing result;Vec direction=sub(end,start);
  for(size_t i=0;i<ports.size();i++){
   if(dot(direction,ports[i].normal)<=0)continue;
   for(auto&t:ports[i].tri){Vec e1=sub(t.b,t.a),e2=sub(t.c,t.a),h=cross(direction,e2);double det=dot(e1,h);if(std::abs(det)<1e-38)continue;
    double inv=1/det;Vec s=sub(start,t.a);double u=inv*dot(s,h);if(u< -1e-10||u>1+1e-10)continue;Vec q=cross(s,e1);double v=inv*dot(direction,q);if(v< -1e-10||u+v>1+1e-10)continue;double f=inv*dot(e2,q);
    if(f>0&&f<=1+1e-10&&f<result.fraction){result={int(i),std::min(f,1.),add(start,mul(direction,f))};}
   }
  }return result;
 }
};
class PositionSampler {
 const frozen::FlowField& field;const Geometry& geometry;std::vector<Tri>tri;std::vector<double>cdf;double total=0,max_un=0;Vec inward;double clearance,pair_clearance;int attempts;bool weighted;
 static double uniform(std::mt19937_64&r){return std::generate_canonical<double,53>(r);}
 public:
 PositionSampler(std::string root,const frozen::FlowField&f,const Geometry&g,double c,double pc,int tries,bool w):field(f),geometry(g),clearance(c),pair_clearance(pc),attempts(tries),weighted(w){
  inward=mul(g.ports[0].normal,-1);std::ifstream in(root+"/geometry/injection_triangles.txt");Tri t;
  while(in>>t.a[0]>>t.a[1]>>t.a[2]>>t.b[0]>>t.b[1]>>t.b[2]>>t.c[0]>>t.c[1]>>t.c[2]){tri.push_back(t);total+=.5*norm(cross(sub(t.b,t.a),sub(t.c,t.a)));cdf.push_back(total);}
  if(tri.empty())throw std::runtime_error("EMPTY_INJECTION_PLANE");
  for(size_t i=0;i<f.velocity.size();i+=3)max_un=std::max(max_un,dot({f.velocity[i],f.velocity[i+1],f.velocity[i+2]},inward));
  if(!(max_un>0))throw std::runtime_error("NO_POSITIVE_INWARD_FLOW");
 }
 bool place(Candidate const&c,Vec&x,std::mt19937_64&r,const std::vector<Admission>&newly,const std::vector<rigid::Particle>&active)const{
  int geometrically_admissible=0;
  for(int n=0;n<attempts;n++){
   size_t i=std::lower_bound(cdf.begin(),cdf.end(),uniform(r)*total)-cdf.begin();auto&t=tri.at(i);double s=std::sqrt(uniform(r)),v=uniform(r);
   x=add(mul(t.a,1-s),add(mul(t.b,s*(1-v)),mul(t.c,s*v)));auto q=frozen::FlowFieldSampler(field).query(x);
   if(q.status!=frozen::Status::VALID)continue;double un=dot(q.u,inward);if(un<=0)continue;
   // Convex trilinear interpolation is bounded by the maximum native-node u.n.
   if(weighted&&uniform(r)*max_un>un)continue;
   if(geometry.wall.distance(x)<c.radius+clearance)continue;
   geometrically_admissible++;
   bool valid=true;for(auto&p:active)if(norm(sub(x,p.x))<c.radius+p.a+pair_clearance){valid=false;break;}
   if(valid)for(auto&p:newly)if(norm(sub(x,p.x))<c.radius+p.candidate.radius+pair_clearance){valid=false;break;}
   if(valid)return true;
  }
  // A finite random search cannot prove size inadmissibility or congestion if
  // it never even tested a geometrically admissible point. Do not mislabel it.
  if(!geometrically_admissible)throw std::runtime_error("POSITION_SEARCH_BUDGET_EXHAUSTED");
  return false;
 }
};
}

#include "surface_patch_topology.hpp"
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iomanip>
namespace workflow {
using namespace rigid;
static Vec cross3(Vec a,Vec b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
static std::pair<int,int> key(int a,int b){return std::minmax(a,b);}
SurfacePatchTopology::SurfacePatchTopology(const std::string&path,const std::string&region_path){
 std::ifstream in(path,std::ios::binary);char h[80];uint32_t count=0;in.read(h,80);in.read(reinterpret_cast<char*>(&count),4);
 if(!in||!count||count>10000000)throw std::runtime_error("STOP_WALL_SURFACE_TOPOLOGY_INVALID");
 std::vector<std::array<Vec,3>> raw(count);std::set<Vec>unique;
 for(auto&t:raw){float a[12];uint16_t attr;in.read(reinterpret_cast<char*>(a),48);in.read(reinterpret_cast<char*>(&attr),2);if(!in)throw std::runtime_error("TRUNCATED_WALL_STL");for(int k=0;k<3;k++){t[k]={a[3+3*k],a[4+3*k],a[5+3*k]};for(double x:t[k])if(!std::isfinite(x))throw std::runtime_error("NONFINITE_WALL_VERTEX");unique.insert(t[k]);}}
 vertices.assign(unique.begin(),unique.end());std::map<Vec,int>ids;for(size_t i=0;i<vertices.size();i++)ids[vertices[i]]=i;
 std::ifstream reg(region_path);faces.resize(count);
 for(size_t i=0;i<count;i++){auto&f=faces[i];f.original_id=i;if(!region_path.empty()&&!(reg>>f.region))throw std::runtime_error("WALL_REGION_MAP_MISMATCH");for(int k=0;k<3;k++)f.vertices[k]=ids.at(raw[i][k]);f.center=mul(add(add(raw[i][0],raw[i][1]),raw[i][2]),1./3);Vec c=cross3(sub(raw[i][1],raw[i][0]),sub(raw[i][2],raw[i][0]));double l=norm(c);if(!(l>0))throw std::runtime_error("STOP_WALL_SURFACE_TOPOLOGY_INVALID");f.area=l/2;f.normal=mul(c,1/l);}
 auto canonical=[](const SurfaceFace&f){auto a=f.vertices;std::sort(a.begin(),a.end());return a;};
 std::sort(faces.begin(),faces.end(),[&](auto&a,auto&b){return canonical(a)<canonical(b);});
 incident.resize(vertices.size());adjacency.resize(count);
 for(size_t i=0;i<count;i++){if(i&&canonical(faces[i])==canonical(faces[i-1]))throw std::runtime_error("DUPLICATE_WALL_TRIANGLE");auto f=faces[i].vertices;for(int k=0;k<3;k++){incident[f[k]].push_back(i);edges[key(f[k],f[(k+1)%3])].incident.push_back(i);}}
 for(auto&entry:edges){auto&e=entry.second;auto const&inc=e.incident;if(inc.size()==1)boundary_edges++;if(inc.size()>2){nonmanifold_edges++;bad_faces.insert(inc.begin(),inc.end());}if(inc.size()!=2)continue;
  int i=inc[0],j=inc[1];e.dihedral=std::acos(std::clamp(dot(faces[i].normal,faces[j].normal),-1.,1.))*180/std::acos(-1.);int orientation=0;
  for(int f:inc)for(int k=0;k<3;k++){int a=faces[f].vertices[k],b=faces[f].vertices[(k+1)%3];if(key(a,b)==entry.first)orientation+=(a<b?1:-1);}
  e.consistent=orientation==0;if(!e.consistent)orientation_conflicts++;
  adjacency[i].push_back(j);adjacency[j].push_back(i);
 }
 order.resize(count);std::iota(order.begin(),order.end(),0);build(0,count);
}
int SurfacePatchTopology::build(int b,int e){Node n;n.begin=b;n.end=e;n.lo={INFINITY,INFINITY,INFINITY};n.hi={-INFINITY,-INFINITY,-INFINITY};for(int i=b;i<e;i++)for(int v:faces[order[i]].vertices)for(int k=0;k<3;k++){n.lo[k]=std::min(n.lo[k],vertices[v][k]);n.hi[k]=std::max(n.hi[k],vertices[v][k]);}int id=nodes.size();nodes.push_back(n);if(e-b>8){int axis=0;for(int k=1;k<3;k++)if(n.hi[k]-n.lo[k]>n.hi[axis]-n.lo[axis])axis=k;int mid=(b+e)/2;std::nth_element(order.begin()+b,order.begin()+mid,order.begin()+e,[&](int i,int j){return faces[i].center[axis]<faces[j].center[axis];});int l=build(b,mid),r=build(mid,e);nodes[id].left=l;nodes[id].right=r;}return id;}
static double box_distance(Vec lo,Vec hi,Vec x){double d=0;for(int k=0;k<3;k++){double z=std::max({lo[k]-x[k],x[k]-hi[k],0.});d+=z*z;}return std::sqrt(d);}
SurfaceCandidate SurfacePatchTopology::closest(int id,Vec x)const{
 auto const&f=faces[id];Vec a=vertices[f.vertices[0]],b=vertices[f.vertices[1]],c=vertices[f.vertices[2]],ab=sub(b,a),ac=sub(c,a),ap=sub(x,a);double aa=dot(ab,ab),bb=dot(ab,ac),cc=dot(ac,ac),pa=dot(ap,ab),pc=dot(ap,ac),det=aa*cc-bb*bb;
 double v=(cc*pa-bb*pc)/det,w=(aa*pc-bb*pa)/det;SurfaceCandidate q;q.id=id;q.original_id=f.original_id;q.region=f.region;q.distance=INFINITY;
 if(v>=0&&w>=0&&v+w<=1){q.point=add(add(a,mul(ab,v)),mul(ac,w));q.barycentric={1-v-w,v,w};q.distance=norm(sub(x,q.point));}
 for(int k=0;k<3;k++){Vec u=vertices[f.vertices[k]],d=sub(vertices[f.vertices[(k+1)%3]],u);double t=std::clamp(dot(sub(x,u),d)/dot(d,d),0.,1.);Vec p=add(u,mul(d,t));double dist=norm(sub(x,p));if(dist<q.distance){q.distance=dist;q.point=p;q.barycentric={0,0,0};q.barycentric[k]=1-t;q.barycentric[(k+1)%3]=t;}}
 return q;
}
void SurfacePatchTopology::nearest(int id,Vec x,double&d)const{auto const&n=nodes[id];if(box_distance(n.lo,n.hi,x)>d)return;if(n.left<0){for(int i=n.begin;i<n.end;i++)d=std::min(d,closest(order[i],x).distance);return;}int a=n.left,b=n.right;if(box_distance(nodes[b].lo,nodes[b].hi,x)<box_distance(nodes[a].lo,nodes[a].hi,x))std::swap(a,b);nearest(a,x,d);nearest(b,x,d);}
void SurfacePatchTopology::collect(int id,Vec x,double d,std::vector<SurfaceCandidate>&out)const{auto const&n=nodes[id];if(box_distance(n.lo,n.hi,x)>d)return;if(n.left<0){for(int i=n.begin;i<n.end;i++){auto q=closest(order[i],x);if(q.distance<=d)out.push_back(q);}return;}collect(n.left,x,d,out);collect(n.right,x,d,out);}
double SurfacePatchTopology::dihedral(int i,int j)const{auto a=faces[i].vertices;for(int k=0;k<3;k++){auto const&e=edges.at(key(a[k],a[(k+1)%3]));if(e.incident.size()==2&&e.consistent&&std::find(e.incident.begin(),e.incident.end(),j)!=e.incident.end())return e.dihedral;}return INFINITY;}
std::set<int> SurfacePatchTopology::one_ring(int face,int vertex,double theta)const{std::set<int>found{face};std::vector<int>stack{face};while(!stack.empty()){int i=stack.back();stack.pop_back();for(int j:adjacency[i])if(!found.count(j)&&std::find(incident[vertex].begin(),incident[vertex].end(),j)!=incident[vertex].end()&&dihedral(i,j)<=theta){found.insert(j);stack.push_back(j);}}return found;}
SurfaceQuery SurfacePatchTopology::query(Vec x,double factor,double theta)const{
 if(!(factor>0&&theta>0&&theta<90))throw std::runtime_error("INVALID_SURFACE_NUMERICAL_PARAMETERS");
 SurfaceQuery q;q.position=x;q.minimum_distance=INFINITY;nearest(0,x,q.minimum_distance);q.tie_tolerance=factor*64*std::numeric_limits<double>::epsilon()*std::max(norm(x),1e-6);collect(0,x,q.minimum_distance+q.tie_tolerance,q.candidates);std::sort(q.candidates.begin(),q.candidates.end(),[](auto&a,auto&b){return a.id<b.id;});
 bool nonmanifold=false;
 for(auto&c:q.candidates){nonmanifold|=bad_faces.count(c.id);std::vector<int>zero;for(int k=0;k<3;k++)if(c.barycentric[k]<=1e-9)zero.push_back(k);std::set<int>support{c.id};
  if(zero.size()>=2){c.feature="VERTEX";int k=std::max_element(c.barycentric.begin(),c.barycentric.end())-c.barycentric.begin();c.vertex=faces[c.id].vertices[k];support=one_ring(c.id,c.vertex,theta);}
  else if(zero.size()==1){c.feature="EDGE";for(int k=0;k<3;k++)if(k!=zero[0])c.edge.push_back(faces[c.id].vertices[k]);std::sort(c.edge.begin(),c.edge.end());auto const&e=edges.at(key(c.edge[0],c.edge[1]));for(int j:e.incident)if(j!=c.id&&dihedral(c.id,j)<=theta)support.insert(j);}
  c.support.assign(support.begin(),support.end());
 }
 if(nonmanifold){q.status="FAIL_NONMANIFOLD";return q;}
 std::vector<int>parent(q.candidates.size());std::iota(parent.begin(),parent.end(),0);auto root=[&](int i){while(parent[i]!=i)i=parent[i];return i;};
 for(size_t i=0;i<parent.size();i++)for(size_t j=i+1;j<parent.size();j++){auto&a=q.candidates[i];auto&b=q.candidates[j];bool same=std::binary_search(a.support.begin(),a.support.end(),b.id)||std::binary_search(b.support.begin(),b.support.end(),a.id);if(same||dihedral(a.id,b.id)<=theta)parent[root(j)]=root(i);}
 std::map<int,std::vector<int>>groups;for(size_t i=0;i<parent.size();i++)groups[root(i)].push_back(i);
 for(auto const&entry:groups){SurfaceCluster c;std::set<int>support,vertex_ids;bool vertex_only=true,edge=false;double pweight=0;
  for(int i:entry.second){auto const&a=q.candidates[i];c.candidates.push_back(a.id);support.insert(a.support.begin(),a.support.end());vertex_only&=a.feature=="VERTEX";edge|=a.feature=="EDGE";if(a.vertex>=0)vertex_ids.insert(a.vertex);double w=faces[a.id].area;pweight+=w;c.point=add(c.point,mul(a.point,w));}
  c.point=mul(c.point,1/pweight);c.support.assign(support.begin(),support.end());c.id=*support.begin();int vertex=vertex_only&&vertex_ids.size()==1?*vertex_ids.begin():-1;double total=0;
  for(int i:support){auto const&f=faces[i];double w=f.area;if(vertex>=0){int k=std::find(f.vertices.begin(),f.vertices.end(),vertex)-f.vertices.begin();Vec a=sub(vertices[f.vertices[(k+1)%3]],vertices[vertex]),b=sub(vertices[f.vertices[(k+2)%3]],vertices[vertex]);w=std::acos(std::clamp(dot(a,b)/(norm(a)*norm(b)),-1.,1.));}c.normal=add(c.normal,mul(f.normal,w));total+=w;}
  double length=norm(c.normal);if(length<1e-14*total){q.status="STOP_WALL_NORMAL_ORIENTATION_INVALID";return q;}c.normal=mul(c.normal,1/length);double orientation=dot(c.normal,sub(x,c.point));if(orientation<0)c.normal=mul(c.normal,-1);if(!(std::abs(orientation)>0)){q.status="STOP_WALL_NORMAL_ORIENTATION_INVALID";return q;}
  c.feature=vertex>=0?"SMOOTH_VERTEX":(edge?"SMOOTH_EDGE":(entry.second.size()==1?"SINGLE_FACE":"SINGLE_SURFACE_PATCH"));q.clusters.push_back(c);
 }
 std::sort(q.clusters.begin(),q.clusters.end(),[](auto&a,auto&b){return a.id<b.id;});q.status=q.clusters.size()>1?"TRUE_MULTI_SURFACE":q.clusters.at(0).feature;
 std::set<int>support;for(auto&c:q.clusters)support.insert(c.support.begin(),c.support.end());q.min_dihedral=INFINITY;for(int i:support)for(int j:adjacency[i])if(support.count(j)){double a=dihedral(i,j);q.min_dihedral=std::min(q.min_dihedral,a);q.max_dihedral=std::max(q.max_dihedral,a);}if(!std::isfinite(q.min_dihedral))q.min_dihedral=0;return q;
}
SurfaceProjection project_surface_velocity(Vec raw,const std::vector<Vec>&N){
 SurfaceProjection best;best.objective=INFINITY;double tol=2e-12*std::max(norm(raw),1e-30);
 auto evaluate=[&](std::vector<int>active){int n=active.size();double a[3][4]{};for(int i=0;i<n;i++){for(int j=0;j<n;j++)a[i][j]=dot(N[active[i]],N[active[j]]);a[i][n]=-dot(N[active[i]],raw);}for(int k=0;k<n;k++){int p=k;for(int i=k+1;i<n;i++)if(std::abs(a[i][k])>std::abs(a[p][k]))p=i;if(std::abs(a[p][k])<1e-12)return;for(int j=k;j<=n;j++)std::swap(a[p][j],a[k][j]);double d=a[k][k];for(int j=k;j<=n;j++)a[k][j]/=d;for(int i=0;i<n;i++)if(i!=k){double v=a[i][k];for(int j=k;j<=n;j++)a[i][j]-=v*a[k][j];}}
  Vec used=raw;std::vector<double>lambda;for(int i=0;i<n;i++){if(a[i][n]<-tol)return;lambda.push_back(a[i][n]);used=add(used,mul(N[active[i]],a[i][n]));}for(Vec normal:N)if(dot(normal,used)<-tol)return;double objective=.5*dot(sub(used,raw),sub(used,raw));if(objective<best.objective){best={used,objective,active,lambda};}
 };
 evaluate({});for(size_t i=0;i<N.size();i++){evaluate({int(i)});for(size_t j=i+1;j<N.size();j++){evaluate({int(i),int(j)});for(size_t k=j+1;k<N.size();k++)evaluate({int(i),int(j),int(k)});}}
 if(!std::isfinite(best.objective))throw std::runtime_error("STOP_MULTINORMAL_PROJECTION_FAILED");return best;
}
void write_vec(std::ostream&o,Vec v){o<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';}
static void ints(std::ostream&o,const std::vector<int>&a){o<<'[';for(size_t i=0;i<a.size();i++){if(i)o<<',';o<<a[i];}o<<']';}
void SurfacePatchTopology::write_json(std::ostream&o,const SurfaceQuery&q)const{
 o<<std::setprecision(17)<<"{\"status\":\""<<q.status<<"\",\"position\":";write_vec(o,q.position);o<<",\"minimum_distance\":"<<q.minimum_distance<<",\"tie_tolerance\":"<<q.tie_tolerance<<",\"min_dihedral\":"<<q.min_dihedral<<",\"max_dihedral\":"<<q.max_dihedral<<",\"candidates\":[";
 for(size_t i=0;i<q.candidates.size();i++){if(i)o<<',';auto const&c=q.candidates[i];o<<"{\"id\":"<<c.id<<",\"original_id\":"<<c.original_id<<",\"region\":"<<c.region<<",\"distance\":"<<c.distance<<",\"point\":";write_vec(o,c.point);o<<",\"barycentric\":";write_vec(o,c.barycentric);o<<",\"feature\":\""<<c.feature<<"\",\"vertex\":"<<c.vertex<<",\"edge\":";ints(o,c.edge);o<<",\"support\":";ints(o,c.support);o<<",\"face_normal\":";write_vec(o,faces[c.id].normal);o<<",\"triangle_coordinates\":[";for(int k=0;k<3;k++){if(k)o<<',';write_vec(o,vertices[faces[c.id].vertices[k]]);}o<<"],\"adjacency\":[";for(size_t k=0;k<adjacency[c.id].size();k++){if(k)o<<',';int j=adjacency[c.id][k];double d=dihedral(c.id,j);o<<"{\"id\":"<<j<<",\"region\":"<<faces[j].region<<",\"dihedral\":"<<(std::isfinite(d)?d:180.)<<'}';}o<<"]}";}
 o<<"],\"clusters\":[";for(size_t i=0;i<q.clusters.size();i++){if(i)o<<',';auto const&c=q.clusters[i];o<<"{\"id\":"<<c.id<<",\"feature\":\""<<c.feature<<"\",\"candidate_ids\":";ints(o,c.candidates);o<<",\"support\":";ints(o,c.support);o<<",\"normal\":";write_vec(o,c.normal);o<<",\"point\":";write_vec(o,c.point);o<<'}';}o<<"]}";
}
}

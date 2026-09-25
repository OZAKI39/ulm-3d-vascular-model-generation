#pragma once
#include "rigid_math.hpp"
#include <map>
#include <set>
#include <ostream>
namespace workflow {
using rigid::Vec;
struct SurfaceCandidate {
 int id=-1,original_id=-1,region=-1,vertex=-1;
 double distance=0;
 Vec point{},barycentric{};
 std::string feature="FACE";
 std::vector<int> edge,support;
};
struct SurfaceCluster {
 int id=-1; Vec normal{},point{};
 std::vector<int> candidates,support;
 std::string feature;
};
struct SurfaceQuery {
 Vec position{};double minimum_distance=0,tie_tolerance=0,min_dihedral=0,max_dihedral=0;
 std::string status;
 std::vector<SurfaceCandidate> candidates;
 std::vector<SurfaceCluster> clusters;
};
struct SurfaceFace {std::array<int,3> vertices;Vec normal{},center{};double area=0;int original_id=-1,region=0;};
struct SurfaceEdge {std::vector<int> incident;double dihedral=0;bool consistent=true;};
struct SurfaceProjection {Vec used{};double objective=0;std::vector<int> active;std::vector<double> multipliers;};
class SurfacePatchTopology {
 public:
 explicit SurfacePatchTopology(const std::string& binary_stl,const std::string& regions="");
 SurfaceQuery query(Vec x,double tie_factor=1,double smooth_deg=15)const;
 void write_json(std::ostream&,const SurfaceQuery&)const;
 std::vector<Vec> vertices;
 std::vector<SurfaceFace> faces;
 std::map<std::pair<int,int>,SurfaceEdge> edges;
 std::vector<std::vector<int>> incident,adjacency;
 int nonmanifold_edges=0,orientation_conflicts=0,boundary_edges=0;
 private:
 struct Node {Vec lo,hi;int left=-1,right=-1,begin=0,end=0;};
 std::vector<Node>nodes;std::vector<int>order;std::set<int>bad_faces;
 int build(int,int);
 SurfaceCandidate closest(int,Vec)const;
 void nearest(int,Vec,double&)const;
 void collect(int,Vec,double,std::vector<SurfaceCandidate>&)const;
 double dihedral(int,int)const;
 std::set<int> one_ring(int,int,double)const;
};
SurfaceProjection project_surface_velocity(Vec raw,const std::vector<Vec>& normals);
void write_vec(std::ostream&,Vec);
}

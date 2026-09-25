#pragma once
#include "rigid_math.hpp"
#include "wall_distance.hpp"
#include <memory>
namespace wallv0 {
using rigid::Vec;using Matrix6=std::array<double,36>;using Frame=std::array<double,9>;
Frame frame(Vec normal,const Frame* previous=nullptr);
Matrix6 rotate(const Matrix6& local,const Frame& Q);
class WallResistanceLookup {
 std::vector<double> epsilon,logs;std::vector<Matrix6> excess;
 public:
 explicit WallResistanceLookup(const std::string&h5);
 Matrix6 scaled_excess(double e)const;
 Matrix6 global_excess(double a,double mu,double gap,const Frame&Q,bool enabled=true)const;
 int validity_class(double e)const{return e<.001?0:(e<=.2?1:(e<=5?2:(e<=20?3:4)));}
};
struct Entry {passive::WallQuery geometry;Frame Q;Matrix6 excess;bool active=false;int validity=0;};
class STLNearestWallQuery {
 bool flat;std::unique_ptr<passive::WallDistance> mesh;
 public:
 explicit STLNearestWallQuery(const std::string&path):flat(path=="FLAT") {if(!flat)mesh=std::make_unique<passive::WallDistance>(path);}
 passive::WallQuery query(Vec x,double radius)const;
 bool swept_safe(Vec x,Vec y,double radius)const;
};
std::vector<rigid::WallBlock> blocks(const std::vector<Entry>&stage,const std::vector<Entry>&base,bool hard);
}
extern "C" void* wall_lookup_open(const char*path);
extern "C" void wall_lookup_close(void*handle);
extern "C" int wall_lookup_evaluate(void*handle,int n,const double*epsilon,double*out);
extern "C" int wall_rotation_audit(int n,const double*normals,const double*matrices,const double*velocities,double*frames,double*forces_matrix,double*forces_two_rotations);
extern "C" void* wall_geometry_open(const char*path);
extern "C" void wall_geometry_close(void*handle);
extern "C" int wall_geometry_query(void*handle,int n,const double*positions,const double*radii,double*out);
extern "C" int wall_geometry_swept(void*handle,int n,const double*x,const double*y,const double*radii,int*out);
extern "C" int wall_single_response(void*handle,int n,const double*state,double*velocity);

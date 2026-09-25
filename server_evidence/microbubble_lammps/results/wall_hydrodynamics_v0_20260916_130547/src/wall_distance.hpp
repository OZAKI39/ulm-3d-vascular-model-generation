#pragma once
#include "frozen_flow.hpp"
#include <vector>
namespace passive {
using frozen::Vec;
struct Triangle { Vec a,b,c,center,normal; int original_id=-1; };
struct Node { Vec lo,hi; int left=-1,right=-1,begin=0,end=0; };
struct WallQuery {
 Vec closest{},normal{},triangle_normal{},triangle_normal_fluid{},plane_normal{};
 double distance=0,gap=0,rms_over_a=0,normal_spread_deg=0,nearest_normal_angle_deg=0;
 int triangle_id=-1,patch_triangles=0;bool valid=false,planar=false,ambiguous=false;
};
class WallDistance {
 public:
 explicit WallDistance(const std::string& path);
 double distance(Vec p) const;
 WallQuery query(Vec p,double radius) const;
 bool segment_safe(Vec a,Vec b,double minimum,int depth=0) const;
 private:
 std::vector<Triangle> triangles; std::vector<Node> nodes;
 int build(int begin,int end); void search(int index,Vec p,double& best) const;
 void nearest(int index,Vec p,double&best,int&which,Vec&point) const;
 void patch(int index,Vec p,double r2,std::vector<int>&out) const;
};
double norm(Vec a); Vec minus(Vec a,Vec b); Vec plus_scaled(Vec a,Vec b,double scale);
}

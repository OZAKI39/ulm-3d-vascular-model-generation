#pragma once
#include "frozen_flow.hpp"
#include <vector>
namespace passive {
using frozen::Vec;
struct Triangle { Vec a,b,c,center; };
struct Node { Vec lo,hi; int left=-1,right=-1,begin=0,end=0; };
struct ClosestWallPoint { Vec point{}; double distance=0; bool ambiguous=false; };
class WallDistance {
 public:
 explicit WallDistance(const std::string& path);
 double distance(Vec p) const;
 ClosestWallPoint closest(Vec p) const;
 bool segment_safe_exact(Vec a,Vec b,double minimum) const;
 bool segment_safe(Vec a,Vec b,double minimum,int depth=0) const;
 private:
 std::vector<Triangle> triangles; std::vector<Node> nodes;
 void closest_search(int index,Vec p,ClosestWallPoint& best) const;
 bool segment_search(int index,Vec a,Vec b,double minimum) const;
 int build(int begin,int end); void search(int index,Vec p,double& best) const;
};
double norm(Vec a); Vec minus(Vec a,Vec b); Vec plus_scaled(Vec a,Vec b,double scale);
}

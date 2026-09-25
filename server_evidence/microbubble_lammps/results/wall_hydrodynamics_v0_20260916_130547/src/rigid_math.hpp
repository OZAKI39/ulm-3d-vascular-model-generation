#pragma once
#include "frozen_flow.hpp"
#include <vector>
#include <utility>
namespace rigid {
using frozen::Vec;
using Mat=std::array<double,9>;
using Array=std::vector<double>;
struct Background{Vec u{},omega{},du{},dw{};Mat grad{},strain{};};
class FlowGradientSampler{const frozen::FlowField&field;public:explicit FlowGradientSampler(const frozen::FlowField&f):field(f){}; Background query(Vec x)const;};
struct Particle{long long id;int owner;Vec x;double a;};
struct Coeff{double sq,sh,pu,tw,li,lj;};
Coeff coefficients(double a,double b,double h);
struct Event{int i,j;double gap;Vec n;Coeff c;Vec slip,normal,transverse,force,ti,tj,fn,fs,tpu;bool constraint=false;};
struct System{Array R,rhs;std::vector<Event> events;};
System assemble(const std::vector<Particle>&p,const std::vector<Background>&bg,const std::vector<std::pair<int,int>>&pairs);
struct WallBlock{int particle_index;long long triangle_id;std::array<double,36> excess;Vec base_normal;double base_gap;bool hard;};
struct Solution{Array q;std::vector<Event>events;double residual=0;int iterations=0;std::vector<int> constraints,wall_constraints;};
Solution solve(const std::vector<Particle>&p,const std::vector<Background>&bg,const std::vector<std::pair<int,int>>&pairs,const std::vector<Particle>*base=nullptr,double dt=0,const std::vector<WallBlock>*walls=nullptr);
Vec add(Vec a,Vec b);Vec sub(Vec a,Vec b);Vec mul(Vec a,double s);double dot(Vec a,Vec b);double norm(Vec a);Vec cross(Vec a,Vec b);Vec mv(Mat a,Vec x);
}
extern "C" void rigid_coefficients(int n,const double*in,double*out);
extern "C" int rigid_audit_system(int n,const double*state,int np,const int*pairs,double*matrix,double*rhs,double*solution,double*stats);

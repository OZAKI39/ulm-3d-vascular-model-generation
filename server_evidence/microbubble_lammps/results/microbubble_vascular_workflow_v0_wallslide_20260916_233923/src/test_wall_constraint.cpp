#include "kinematic_wall_constraint.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace rigid;
using namespace workflow;
int main(int argc,char**argv){
 double a=1e-6,m=1e-10,V=1e-4;Vec x{0,0,a+m};
 WallQuery plane=[](Vec x){return passive::ClosestWallPoint{{x[0],x[1],0},x[2],false};};
 WallSegment ps=[](Vec a,Vec b,double d){return std::min(a[2],b[2])>=d;};
 for(Vec raw: {Vec{0,0,-V},Vec{V,0,0},Vec{V,0,-V},Vec{V,0,V}}){auto r=constrain_wall({x,x,raw,a,m,1e-4,false},plane,ps);Vec expected=raw;expected[2]=std::max(raw[2],0.);assert(norm(sub(expected,r.constrained_velocity))/V<1e-12);assert(r.safe);}
 // Rotate the complete geometry and velocity by a fixed orthonormal rotation.
 Vec n=mul(Vec{1,2,3},1/std::sqrt(14.));Vec t=mul(Vec{-2,1,0},1/std::sqrt(5.));Vec xr=mul(n,a+m),vr=mul(sub(t,n),V);
 WallQuery rotated=[&](Vec x){double d=dot(x,n);return passive::ClosestWallPoint{sub(x,mul(n,d)),d,false};};
 WallSegment rs=[&](Vec x,Vec y,double d){return std::min(dot(x,n),dot(y,n))>=d;};
 auto r=constrain_wall({xr,xr,vr,a,m,1e-4,false},rotated,rs);assert(norm(sub(r.constrained_velocity,mul(t,V)))/V<1e-12);assert(std::abs(dot(r.constrained_velocity,n))<1e-16);
 if(argc>1){passive::WallDistance wall(argv[1]);assert(wall.segment_safe_exact({0,0,a+m},{V*1e-4,0,a+m},a+m));assert(!wall.segment_safe_exact({0,0,a+m},{0,0,a},a+m));}
 std::ofstream csv("SYNTHETIC_CURVED.csv");csv<<std::setprecision(17)<<"shape,dt_max,step,time,x,y,z,gap,path_length,constraint_events,dt,max_correction,total_correction,retries\n";
 for(std::string shape:{"cylinder","sphere"})for(double dtmax:{1e-4,5e-5,2.5e-5}){
  double R=5e-6,rc=R-a-m;Vec x{rc,0,0};double time=0,path=0,total=0,maxcor=0;int step=0,events=0,retries=0;
  auto rad=[&](Vec x){return shape=="cylinder"?std::hypot(x[0],x[1]):norm(x);};
  WallQuery query=[&](Vec x){double r=rad(x);Vec p=mul(x,R/r);if(shape=="cylinder")p[2]=x[2];return passive::ClosestWallPoint{p,R-r,false};};
  WallSegment segment=[&](Vec x,Vec y,double d){return std::max(rad(x),rad(y))<=R-d;};
  auto velocity=[&](Vec x){Vec radial=mul(x,1/rad(x));if(shape=="cylinder")radial[2]=0;Vec tangent{-radial[1],radial[0],0};return mul(add(radial,tangent),V);};
  while(time<.01-1e-14){double dt=std::min(dtmax,.01-time);bool ok=false;WallConstraintResult f,g;
   for(int attempt=0;attempt<24;attempt++){f=constrain_wall({x,x,velocity(x),a,m,dt/2,attempt>0},query,segment);if(f.safe){g=constrain_wall({f.endpoint,x,velocity(f.endpoint),a,m,dt,attempt>0},query,segment);if(g.safe){ok=true;break;}}dt/=2;retries++;}
   assert(ok);path+=norm(sub(g.endpoint,x));x=g.endpoint;time+=dt;step++;events+=f.active+g.active;total+=f.position_projection_distance+g.position_projection_distance;maxcor=std::max({maxcor,f.position_projection_distance,g.position_projection_distance});assert(R-rad(x)-a>=m-2e-14);assert(maxcor<=m/4);
   csv<<shape<<','<<dtmax<<','<<step<<','<<time;for(double v:x)csv<<','<<v;csv<<','<<R-rad(x)-a<<','<<path<<','<<events<<','<<dt<<','<<maxcor<<','<<total<<','<<retries<<'\n';
  }
 }
 std::cout<<"PLANE_A_B_C_D_E_PASS; CURVED_F_COMPLETED; REFINEMENT_DATA_WRITTEN\n";
}

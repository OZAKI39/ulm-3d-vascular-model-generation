#include "kinematic_wall_constraint.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
using namespace workflow;using namespace rigid;
int main(int argc,char**argv){try{if(argc!=2)return 2;std::string root=argv[1];std::ofstream out(root+"/raw/REMOTE_SYNTHETIC_TRAJECTORIES.csv");out<<std::setprecision(17)<<"mesh,dt_max,time,x,y,z,dt,stage,swept_gap,correction,candidates,clusters\n";
 for(std::string name:{"plane0","plane1","plane2","plane0_reordered","cylinder","cylinder_reordered","sphere","sphere_reordered"}){
 bool plane=name.find("plane")==0,reorder=name.find("reordered")!=std::string::npos;std::string path=root+"/raw/synthetic/"+name+".stl";SurfacePatchTopology top(path);passive::WallDistance wall(path);
 for(double dtmax:{1e-4,5e-5,2.5e-5,1.25e-5}){if((plane||reorder)&&dtmax!=1e-4)continue;double a=1e-6,margin=1e-10,T=plane?.02:.01,t=0;Vec x;
 if(plane)x={-1e-6,0,a+margin+1e-17};else{Vec d{std::cos(.021),std::sin(.021),0};double lo=0,hi=8e-6;for(int k=0;k<70;k++){double m=(lo+hi)/2;if(wall.distance(mul(d,m))>a+margin+1e-17)lo=m;else hi=m;}x=mul(d,lo);}
 auto raw=[&](Vec p){if(plane)return Vec{1e-4,0,-2e-5};p[2]=0;p=mul(p,1/norm(p));return add(Vec{-1e-4*p[1],1e-4*p[0],0},mul(p,2e-5));};int steps=0,retries=0;std::set<int>faces;
 while(t<T-1e-15){double dt=std::min(dtmax,T-t);bool offset=false,ok=false;WallConstraintResult r0,r1;for(int attempt=0;attempt<40&&dt>=1e-14;attempt++){r0=constrain_surface_wall({x,x,raw(x),a,margin,dt/2,offset},wall,top);if(r0.safe){r1=constrain_surface_wall({r0.endpoint,x,raw(r0.endpoint),a,margin,dt,offset},wall,top);if(r1.safe){ok=true;break;}}dt/=2;offset=true;retries++;}if(!ok)throw std::runtime_error("SYNTHETIC_TIME_REFINEMENT_FAILED "+name);int stage=0;for(auto const&r:{r0,r1}){double gap=wall.segment_distance_exact(x,r.endpoint)-a;if(gap<margin-3e-20||r.surface.clusters.size()!=1||r.position_projection_distance>margin/4)throw std::runtime_error("SYNTHETIC_SAFETY_FAILED");for(auto&c:r.surface.candidates)faces.insert(c.id);out<<name<<','<<dtmax<<','<<t+dt;for(double z:r.endpoint)out<<','<<z;out<<','<<dt<<','<<stage++<<','<<gap<<','<<r.position_projection_distance<<','<<r.surface.candidates.size()<<','<<r.surface.clusters.size()<<'\n';}x=r1.endpoint;t+=dt;if(++steps>100000)throw std::runtime_error("SYNTHETIC_STEP_CAP");}
 std::cout<<std::setprecision(17)<<name<<" dt="<<dtmax<<" time="<<t<<" steps="<<steps<<" retries="<<retries<<" faces="<<faces.size()<<" PASS\n";
 }}return 0;}catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 2;}}

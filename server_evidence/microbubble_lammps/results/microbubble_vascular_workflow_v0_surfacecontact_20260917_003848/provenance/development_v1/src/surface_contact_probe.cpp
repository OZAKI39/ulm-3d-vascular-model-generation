#include "kinematic_wall_constraint.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
using namespace workflow;
int main(int argc,char**argv){try{if(argc<2)return 2;SurfacePatchTopology top(argv[1],argc>2?argv[2]:"");passive::WallDistance wall(argv[1]);std::cout<<std::setprecision(17)<<std::unitbuf;std::string line;
 while(std::getline(std::cin,line)){std::istringstream in(line);char mode;Vec x,v;double factor,theta;in>>mode;for(double&a:x)in>>a;for(double&a:v)in>>a;in>>factor>>theta;if(!in)return 3;
 if(mode=='Q'){auto q=top.query(x,factor,theta);std::cout<<"{\"query\":";top.write_json(std::cout,q);if(!q.clusters.empty()){std::vector<Vec>N;for(auto&c:q.clusters)N.push_back(c.normal);auto p=project_surface_velocity(v,N);std::cout<<",\"used\":";write_vec(std::cout,p.used);std::cout<<",\"objective\":"<<p.objective<<",\"active\":[";for(size_t i=0;i<p.active.size();i++){if(i)std::cout<<',';std::cout<<p.active[i];}std::cout<<']';}std::cout<<"}\n";}
 else if(mode=='C'){double a,margin,dt;int offset;Vec base;in>>a>>margin>>dt>>offset;for(double&b:base)in>>b;auto r=constrain_surface_wall({x,base,v,a,margin,dt,bool(offset)},wall,top,factor,theta);std::cout<<"{\"query\":";top.write_json(std::cout,r.surface);std::cout<<",\"used\":";write_vec(std::cout,r.constrained_velocity);std::cout<<",\"endpoint\":";write_vec(std::cout,r.endpoint);std::cout<<",\"correction\":"<<r.position_projection_distance<<",\"safe\":"<<r.safe<<",\"raw_safe\":"<<r.raw_segment_safe<<",\"mode\":\""<<r.mode<<"\",\"fraction\":"<<r.first_contact_fraction<<",\"swept_gap\":"<<wall.segment_distance_exact(base,r.endpoint)-a<<"}\n";}
 else return 4;
 }return 0;}catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 2;}}

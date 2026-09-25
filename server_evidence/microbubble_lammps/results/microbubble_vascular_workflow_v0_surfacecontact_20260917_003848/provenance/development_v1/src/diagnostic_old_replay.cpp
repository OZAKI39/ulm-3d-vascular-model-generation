// Diagnostic only: links the immutable old wallslide library, not the new contact implementation.
#include "kinematic_wall_constraint.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
using namespace rigid;
void vec(std::ostream&o,Vec x){o<<'['<<x[0]<<','<<x[1]<<','<<x[2]<<']';}
int main(int argc,char**argv){if(argc!=2)return 2;std::string root=argv[1];std::ifstream input(root+"/raw/OLD_LONG_LAST_ACCEPTED.txt");std::vector<Particle>p;Particle t;while(input>>t.id>>t.a>>t.x[0]>>t.x[1]>>t.x[2]){t.owner=0;p.push_back(t);}frozen::FlowField f(root+"/fields/FROZEN_FLOW_FIELD_V0.h5");FlowGradientSampler sampler(f);passive::WallDistance wall(root+"/geometry/WALL_ONLY.stl");std::vector<std::pair<int,int>>pairs;for(int i=0;i<int(p.size());i++)for(int j=i+1;j<int(p.size());j++)pairs.emplace_back(i,j);
 auto bg=[&](const std::vector<Particle>&p){std::vector<Background>b;for(auto&q:p)b.push_back(sampler.query(q.x));return b;};
 std::ofstream log(root+"/raw/OLD_LONG_DIAGNOSTIC_QUERIES.jsonl");log<<std::setprecision(17);double dt=1e-4;bool prior_wall_retry=false;
 for(int attempt=0;attempt<=24&&dt>=1e-10;attempt++){
  auto b=bg(p);auto cur=solve(p,b,pairs);double allowed=dt;for(size_t i=0;i<p.size();i++){double v=norm({cur.q[6*i],cur.q[6*i+1],cur.q[6*i+2]});if(v>0)allowed=std::min(allowed,.25*f.dx/v);}if(allowed<dt*(1-1e-12)){dt=allowed;continue;}
  auto mid=p,end=p;bool safe=true;Vec query_point{},raw{},eval{};passive::ClosestWallPoint closest;int stage=0,index=0;double h=dt/2;
  try{
   for(stage=0;stage<2&&safe;stage++){
    auto const&e=stage==0?p:mid;auto q=solve(e,stage==0?b:bg(mid),pairs,&p,stage==0?dt/2:dt);h=stage==0?dt/2:dt;
    for(index=0;index<int(p.size());index++){
     raw={q.q[6*index],q.q[6*index+1],q.q[6*index+2]};eval=e[index].x;
     auto query=[&](Vec x){query_point=x;closest=wall.closest(x);log<<"{\"attempt\":"<<attempt<<",\"stage\":"<<stage<<",\"particle_id\":"<<p[index].id<<",\"dt\":"<<dt<<",\"query_point\":";vec(log,x);log<<",\"distance\":"<<closest.distance<<",\"closest_point\":";vec(log,closest.point);log<<",\"ambiguous\":"<<(closest.ambiguous?"true":"false")<<"}\n";return closest;};
     auto result=workflow::constrain_wall({eval,p[index].x,raw,p[index].a,1e-10,h,prior_wall_retry},query,[&](Vec a,Vec b,double d){return wall.segment_safe_exact(a,b,d);});
     (stage==0?mid:end)[index].x=result.endpoint;safe=safe&&result.safe;
    }
   }
  }catch(std::exception const&e){
   std::ofstream out(root+"/validation/OLD_LONG_FAILURE_REPLAY.json");out<<std::setprecision(17)<<"{\"status\":\"REPRODUCED\",\"reason\":\""<<e.what()<<"\",\"time\":0.49412276564811575,\"step\":18134,\"stage\":"<<stage<<",\"attempt\":"<<attempt<<",\"dt\":"<<dt<<",\"stage_dt\":"<<h<<",\"particle_id\":"<<p[index].id<<",\"radius\":"<<p[index].a<<",\"evaluation_point\":";vec(out,eval);out<<",\"base\":";vec(out,p[index].x);out<<",\"raw_velocity\":";vec(out,raw);out<<",\"query_point\":";vec(out,query_point);out<<",\"old_closest_point\":";vec(out,closest.point);out<<",\"old_distance\":"<<closest.distance<<",\"old_ambiguous\":"<<(closest.ambiguous?"true":"false")<<"}\n";out.flush();log.flush();std::cout<<"REPRODUCED "<<e.what()<<" stage="<<stage<<" id="<<p[index].id<<" attempt="<<attempt<<'\n';return 0;
  }
  if(safe){std::cerr<<"OLD_FAILURE_NOT_REPRODUCED\n";return 3;}prior_wall_retry=true;dt*=.5;
 }
 std::cerr<<"OLD_FAILURE_NOT_REPRODUCED\n";return 4;
}

#include "adaptive_injection.hpp"
#include <iostream>
#include <cassert>
using namespace workflow;
int main(){
 NumberFluxProvider physical(2.,3.,false),direct(2.,7.,true);assert(physical.number_flux(9)==6&&direct.number_flux(9)==7);
 auto draw=[](long long i){return std::array<double,2>{.01*i,1.+(i%4==0?2.:0.)};};
 auto all=[](Candidate const&){return true;};
 auto place=[](Candidate const&c,std::array<double,3>&x,auto&,auto const&){x={double(c.id),0,0};return true;};
 int tests=0;
 {Controller c;std::vector<InjectionEvent>e;auto admitted=c.advance(0,.5,2,draw,all,place,e);assert(admitted.size()==1&&c.target==1&&c.admitted==1);c.advance(.5,.5,2,draw,all,place,e);assert(c.target==2&&c.admitted==2);tests++;}
 for(bool source:{false,true})for(double lambda:{.3,1.,3.7}){
  Controller c;c.source_basis=source;double t=0,expected=0;
  for(int i=0;i<30;i++){double dt=(i%3+1)*.7;std::vector<InjectionEvent>e;auto a=c.advance(t,dt,lambda,draw,all,place,e);expected+=dt*lambda;t+=dt;assert(c.admitted==std::floor(expected+1e-10));assert(std::abs(c.target-expected)<1e-12);}
  tests++;
 }
 // Time dependent rates and varying timesteps are supplied at the API boundary.
 Controller variable;double expected=0;
 for(int i=0;i<50;i++){std::vector<InjectionEvent>e;double dt=.1+.003*i,rate=2+.1*i;expected+=dt*rate;variable.advance(i,dt,rate,draw,all,place,e);}
 assert(variable.admitted==std::floor(expected+1e-10));tests++;
 for(bool source:{true,false}){
  Controller c;c.source_basis=source;auto feasible=[](Candidate const&c){return c.radius<2;};std::vector<InjectionEvent> e;
  c.advance(0,10,1,draw,feasible,place,e);assert(c.pending.empty());assert(c.size_rejected>0);
  if(source){assert(c.source_drawn==10&&c.admitted==7);}else{assert(c.admitted==10&&c.source_drawn==14);}tests++;
 }
 Controller congested;auto blocked=[](Candidate const&,auto&,auto&,auto const&){return false;};std::vector<InjectionEvent>e;
 congested.advance(0,5,1,draw,all,blocked,e);assert(congested.pending.size()==5);auto first=congested.pending.front();
 congested.advance(5,.1,1,draw,all,blocked,e);assert(congested.pending.front().id==first.id&&congested.pending.front().radius==first.radius);assert(congested.capacity_exceeded(10,4,100));
 auto admitted=congested.advance(5.1,.1,1,draw,all,place,e);assert(admitted.size()==5);for(int i=0;i<5;i++)assert(admitted[i].candidate.id==i+1);tests++;
 // A rejected long trial that draws many particles must not alter the shorter retry.
 Controller base,trial=base;trial.advance(0,10,2,draw,all,place,e);assert(base.source_drawn==0&&base.target==0);
 trial=base;trial.advance(0,.2,2,draw,all,place,e);base=trial;assert(base.source_drawn==0&&base.target==.4);tests++;
 Controller ids;ids.advance(0,8,1,draw,[](Candidate const&c){return c.id%2==0;},place,e);assert(ids.admitted==4&&ids.source_drawn==8);tests++;
 std::cout<<"{\"status\":\"PASS\",\"tests\":"<<tests<<",\"semantics\":[\"SOURCE_POPULATION\",\"ADMITTED_POPULATION\"],\"multiple_per_step\":true,\"FIFO_identity\":true,\"transaction_rollback\":true,\"disclaimer\":\""<<disclaimer<<"\"}\n";
}

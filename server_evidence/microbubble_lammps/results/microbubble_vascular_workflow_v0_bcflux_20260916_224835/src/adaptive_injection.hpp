#pragma once
#include <cmath>
#include <deque>
#include <vector>
#include <string>
#include <array>
#include <stdexcept>
#include <algorithm>
#include <random>
namespace workflow {
constexpr const char* disclaimer="NOT EXPERIMENTAL CONCENTRATION";
class AuthoritativeBoundaryFlowProvider {
 double Q;
 public: explicit AuthoritativeBoundaryFlowProvider(double value):Q(value){if(!std::isfinite(Q)||Q<=0)throw std::runtime_error("INVALID_AUTHORITATIVE_Q");}
 double volume_flux(double time)const{if(!std::isfinite(time))throw std::runtime_error("INVALID_TIME");return Q;}
};
struct NumberFluxProvider {
 AuthoritativeBoundaryFlowProvider flow; double value; bool direct;
 NumberFluxProvider(double Q,double v,bool test):flow(Q),value(v),direct(test){if(!std::isfinite(v)||v<0)throw std::runtime_error("UNSPECIFIED_CONCENTRATION_OR_FLUX");}
 double number_flux(double t)const{return direct?value:value*flow.volume_flux(t);}
};
struct Candidate {long long id;double radius,born,u;};
struct InjectionEvent {long long id;std::string kind;double radius,born,time,u;std::array<double,3> x{};};
struct Admission {Candidate candidate;std::array<double,3> x;};
struct Controller {
 bool source_basis=true;double target=0,compensation=0;long long source_drawn=0,size_rejected=0,admitted=0,id_offset=0;
 std::deque<Candidate> pending;std::mt19937_64 position_rng{20260916};
 // All mutations happen on a trial copy owned by rank zero. Commit only after
 // accepted RK2; rejected trial copies and inserted LAMMPS atoms are discarded.
 template<class Draw,class Feasible,class Place>
 std::vector<Admission> advance(double time,double dt,double lambda,Draw draw,Feasible feasible,Place place,std::vector<InjectionEvent>&events,int draw_budget=1000){
  if(!(dt>0&&lambda>=0&&std::isfinite(lambda)))throw std::runtime_error("INVALID_FLUX_INTERVAL");
  double y=lambda*dt-compensation,z=target+y;compensation=(z-target)-y;target=z;
  std::vector<Admission> out;
  auto attempt=[&](){while(!pending.empty()){
    auto c=pending.front();std::array<double,3>x{};
    if(!place(c,x,position_rng,out)){events.push_back({c.id,"TRANSIENTLY_BLOCKED",c.radius,c.born,time,c.u,{}});break;}
    pending.pop_front();admitted++;out.push_back({c,x});events.push_back({c.id,"ADMITTED",c.radius,c.born,time,c.u,x});
  }};
  attempt(); // existing FIFO candidates first
  long long goal=static_cast<long long>(std::floor(target+1e-10));int draws=0;
  auto covered=[&](){return source_basis?source_drawn:admitted+static_cast<long long>(pending.size());};
  while(covered()<goal){
    if(draws++>=draw_budget)throw std::runtime_error("SOURCE_DRAW_BUDGET_EXCEEDED");
    auto sample=draw(source_drawn);source_drawn++;Candidate c{id_offset+source_drawn,sample[1],time,sample[0]};
    events.push_back({c.id,"SOURCE_DRAWN",c.radius,c.born,time,c.u,{}});
    if(!feasible(c)){size_rejected++;events.push_back({c.id,"SIZE_INADMISSIBLE_AT_INLET",c.radius,c.born,time,c.u,{}});}
    else pending.push_back(c);
  }
  attempt();
  if(source_drawn!=size_rejected+static_cast<long long>(pending.size())+admitted)throw std::runtime_error("FAIL_PARTICLE_ACCOUNTING");
  return out;
 }
 double flux_debt()const{return target-admitted;}
 bool capacity_exceeded(double time,size_t max_pending,double max_age)const{return pending.size()>max_pending||(!pending.empty()&&time-pending.front().born>max_age);}
};
}

#include "flat_wall_engine.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>

int main(int argc,char**argv) {
 try {
  if(argc!=3)throw std::runtime_error("usage: contract TABLE OUTPUT_CSV");
  using namespace rigid;using reducedwall::Status;
  phase2::Plane plane{{0,0,0},{{1,0,0},{0,1,0},{0,0,1}},100};
  phase2::FlatWallEngine wall(argv[1],plane);reducedwall::ReducedWallKernel kernel(argv[1]);
  std::ofstream f(argv[2]);f<<std::setprecision(17)<<"test,value,limit,pass\n";int count=0;
  auto check=[&](const char*name,double value,double limit){bool pass=value<=limit;f<<name<<','<<value<<','<<limit<<','<<pass<<'\n';++count;if(!pass)throw std::runtime_error(name);};
  const double a=std::ldexp(1.,-20);
  for(double e:{.001,.2}) {
   // Direct exact domain endpoints: no coordinate subtraction in this kernel contract probe.
   auto w=kernel.evaluate_shear_rhs({1.,.001,e,plane.frame,{100,0,0},reducedwall::WallShearSource::ANALYTIC});
   check(e==.001?"lower_endpoint":"upper_endpoint",w.status==Status::OK?0:1,0);
  }
  for(double e:{.0009,.2001,21.}) {
   auto w=kernel.evaluate_shear_rhs({1.,.001,e,plane.frame,{100,0,0},reducedwall::WallShearSource::ANALYTIC});
   check("outside_domain_status",std::string(phase2::engine_status(w.status))=="OUTSIDE_CERTIFIED_WALL_DOMAIN"?0:1,0);
  }
  for(double h:{0.,-a})check("invalid_gap",kernel.evaluate_shear_rhs({a,.001,h,plane.frame,{100,0,0},reducedwall::WallShearSource::ANALYTIC}).status==Status::INVALID_GAP?0:1,0);
  for(auto source:{reducedwall::WallShearSource::GRADIENT_PROXY,reducedwall::WallShearSource::PALABOS_VALIDATED})
   check("unsupported_shear_source",kernel.evaluate_shear_rhs({a,.001,.02*a,plane.frame,{100,0,0},source}).status==Status::UNSUPPORTED_SOURCE?0:1,0);
  Vec start{0,0,1.02*a},crossing{a,0,.99*a},touching{a,0,a},safe{a,0,1.02*a};
  const auto original=crossing;
  check("crossing_proposal_rejected",phase2::segment_safe(plane,start,crossing,a)?1:0,0);
  check("contact_proposal_rejected",phase2::segment_safe(plane,start,touching,a)?1:0,0);
  check("safe_proposal_accepted",phase2::segment_safe(plane,start,safe,a)?0:1,0);
  check("proposal_not_repaired",crossing==original?0:1,0);
  check("nan_proposal_rejected",phase2::segment_safe(plane,start,{0,0,std::numeric_limits<double>::quiet_NaN()},a)?1:0,0);
  // Query two genuinely different heights: ensure the adapter evaluates its current argument.
  Particle p{1,0,start,a},m=p;m.x[2]=1.03*a;
  auto wp=wall.evaluate(p),wm=wall.evaluate(m);
  check("requery_changed_gap",wp.resistance.epsilon!=wm.resistance.epsilon?0:1,0);
  check("requery_changed_resistance",wp.resistance.wall_excess_SI!=wm.resistance.wall_excess_SI?0:1,0);
  check("requery_changed_shear_rhs",wp.b_wall_shear!=wm.b_wall_shear?0:1,0);
  std::vector<Particle>ps{p};std::vector<Background>bg{wall.background(p)};
  std::vector<reducedwall::ReducedWallShearContribution>ws{wp};
  auto bulk=assemble(ps,bg,{}),total=assemble(ps,bg,{},&ws);double rhs_error=0,matrix_error=0;
  for(int i=0;i<6;++i){rhs_error=std::max(rhs_error,std::abs(total.rhs[i]-bulk.rhs[i]-wp.b_wall_shear[i])/std::max(std::abs(total.rhs[i]),1e-30));
   for(int j=0;j<6;++j)matrix_error=std::max(matrix_error,std::abs(total.R[6*i+j]-bulk.R[6*i+j]-wp.resistance.wall_excess_SI[6*i+j])/std::max(std::abs(total.R[6*i+j]),1e-30));}
  check("one_wall_rhs_addition",rhs_error,1e-12);check("one_wall_matrix_addition",matrix_error,1e-12);
  auto q=solve(ps,bg,{},nullptr,0,&ws);double err=0;
  for(int i=0;i<6;++i)err=std::max(err,std::abs(q.q[i]-wp.q_target[i])/(i<3?100*a:100));
  check("actual_PCG_recovery",err,1e-9);
  std::cout<<"CONTRACT_PASS checks="<<count<<'\n';return 0;
 }catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}
}

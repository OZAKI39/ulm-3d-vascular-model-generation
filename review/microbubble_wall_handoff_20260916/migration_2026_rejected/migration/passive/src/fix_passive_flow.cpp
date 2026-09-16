#include "fix_passive_flow.hpp"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "error.h"
#include "timer.h"
#include <mpi.h>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace LAMMPS_NS;
using namespace frozen;
using namespace passive;
FixPassiveFlow::FixPassiveFlow(LAMMPS*lmp,int argc,char**argv):Fix(lmp,argc,argv){
 if(argc!=9)error->all(FLERR,"sonovue/passive/flow expects FIELD STRIDE PREFIX STL_OR_NONE MARGIN_M CHECK_PAIRS");
 try{field=std::make_unique<FlowField>(argv[3]);stride=std::stoi(argv[4]);prefix=argv[5];if(std::string(argv[6])!="NONE")wall=std::make_unique<WallDistance>(argv[6]);margin=std::stod(argv[7]);check_pairs=std::stoi(argv[8]);if(stride<1||margin<0)throw std::runtime_error("Invalid passive fix parameters");}catch(std::exception const&e){error->all(FLERR,e.what());}
 nevery=1;time_integrate=1;MPI_Comm_rank(world,&rank);
 trace.open(prefix+"_rank"+std::to_string(rank)+".csv");if(!trace)error->all(FLERR,"Cannot open passive trajectory");trace<<std::setprecision(17)<<"particle_id,time_s,step,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,fluid_ux_m_s,fluid_uy_m_s,fluid_uz_m_s,speed_m_s,fluid_speed_m_s,diameter_um,center_wall_distance_um,surface_wall_gap_um,query_status,owner_rank\n";
}
int FixPassiveFlow::setmask(){return FixConst::INITIAL_INTEGRATE|FixConst::END_OF_STEP|FixConst::POST_RUN;}
bool FixPassiveFlow::sample(Vec x,Vec&u){++queries;auto q=FlowFieldSampler(*field).query(x);if(q.status!=Status::VALID){++invalid;return false;}u=q.u;for(double a:u)if(!std::isfinite(a)){++nonfinite;return false;}return true;}
void FixPassiveFlow::stop(int why){stopped=true;reason=why;timer->force_timeout();}
void FixPassiveFlow::prepare_next(){
 pending.clear();int local_reason=0;
 for(int i=0;i<atom->nlocal;++i)if(atom->mask[i]&groupbit){Proposal p;p.radius=atom->radius[i];for(int a=0;a<3;++a)p.old[a]=atom->x[i][a];Vec k1,k2;
  if(!sample(p.old,k1)){local_reason=std::max(local_reason,3);continue;}
  p.mid=plus_scaled(p.old,k1,.5*update->dt);
  if(wall&&wall->distance(p.mid)<p.radius+margin){local_reason=std::max(local_reason,1);continue;}
  if(!sample(p.mid,k2)){local_reason=std::max(local_reason,3);continue;}
  p.next=plus_scaled(p.old,k2,update->dt);
  if(wall&&!wall->segment_safe(p.old,p.next,p.radius+margin)){local_reason=std::max(local_reason,1);continue;}
  if(!sample(p.next,p.velocity)){local_reason=std::max(local_reason,3);continue;}
  pending.emplace(atom->tag[i],p);
 }
 int global_reason=0;MPI_Allreduce(&local_reason,&global_reason,1,MPI_INT,MPI_MAX,world);
 if(!global_reason&&check_pairs){
  // Only the small real-field population enables this global geometric safety test.
  std::vector<double> send;for(auto const&kv:pending){auto const&p=kv.second;for(double x:p.old)send.push_back(x);for(double x:p.next)send.push_back(x);send.push_back(p.radius);}
  int n=send.size(),ranks;MPI_Comm_size(world,&ranks);std::vector<int>counts(ranks),displs(ranks);MPI_Allgather(&n,1,MPI_INT,counts.data(),1,MPI_INT,world);for(int k=1;k<ranks;++k)displs[k]=displs[k-1]+counts[k-1];std::vector<double>all(displs.back()+counts.back());MPI_Allgatherv(send.data(),n,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,world);
  for(size_t i=0;i<all.size()/7;++i)for(size_t j=i+1;j<all.size()/7;++j){Vec r,v;double vv=0,rv=0;for(int k=0;k<3;++k){r[k]=all[7*i+k]-all[7*j+k];v[k]=(all[7*i+3+k]-all[7*i+k])-(all[7*j+3+k]-all[7*j+k]);vv+=v[k]*v[k];rv+=r[k]*v[k];}double t=vv>0?std::clamp(-rv/vv,0.,1.):0.;if(norm(plus_scaled(r,v,t))<all[7*i+6]+all[7*j+6])global_reason=2;}
 }
 if(global_reason)stop(global_reason);
}
void FixPassiveFlow::setup(int){int bad=0;for(int i=0;i<atom->nlocal;++i)if(atom->mask[i]&groupbit){Vec x,u;for(int a=0;a<3;++a)x[a]=atom->x[i][a];if(!sample(x,u)){bad=3;continue;}for(int a=0;a<3;++a)atom->v[i][a]=u[a];}int global=0;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,world);if(global)stop(global);record(true);if(!stopped&&update->ntimestep<update->laststep)prepare_next();}
void FixPassiveFlow::initial_integrate(int){if(stopped)error->all(FLERR,"Passive step attempted after hard safety stop");for(int i=0;i<atom->nlocal;++i)if(atom->mask[i]&groupbit){auto it=pending.find(atom->tag[i]);if(it==pending.end())error->all(FLERR,"Missing passive proposal");for(int a=0;a<3;++a){atom->x[i][a]=it->second.next[a];atom->v[i][a]=it->second.velocity[a];}}++accepted;}
void FixPassiveFlow::record(bool force){bool write=update->ntimestep!=last_written&&(force||update->ntimestep%stride==0||update->ntimestep==update->laststep);int bad=0;
 for(int i=0;i<atom->nlocal;++i)if(atom->mask[i]&groupbit){Vec x,v,u;for(int a=0;a<3;++a){x[a]=atom->x[i][a];v[a]=atom->v[i][a];if(!std::isfinite(v[a])||!std::isfinite(x[a])){++nonfinite;bad=4;}}
  if(!sample(x,u)){bad=std::max(bad,3);continue;}double wd=wall?wall->distance(x):-1.;if(wall&&wd<atom->radius[i]+margin)bad=std::max(bad,4);
  if(write){double xu[3];domain->unmap(atom->x[i],atom->image[i],xu);trace<<atom->tag[i]<<','<<update->ntimestep*update->dt<<','<<update->ntimestep;for(double a:xu)trace<<','<<a;for(auto const&arr:{v,u})for(double a:arr)trace<<','<<a;trace<<','<<norm(v)<<','<<norm(u)<<','<<2e6*atom->radius[i]<<','<<(wall?wd*1e6:-1.)<<','<<(wall?(wd-atom->radius[i])*1e6:-1.)<<",0,"<<rank<<'\n';}}
 int global=0;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,world);if(global)stop(global);if(write){last_written=update->ntimestep;trace.flush();}}
void FixPassiveFlow::end_of_step(){record();if(!stopped&&update->ntimestep<update->laststep)prepare_next();if(stopped)record(true);}
void FixPassiveFlow::post_run(){record(true);trace.flush();std::ofstream f(prefix+"_state_rank"+std::to_string(rank)+".json");const char*reasons[]={"MAX_PHYSICAL_TIME","WALL_SAFETY_STOP","SPHERE_OVERLAP_SAFETY_STOP","INVALID_QUERY_HARD_STOP","NONFINITE_OR_UNSAFE_STATE"};f<<std::setprecision(17)<<"{\"completed_step\":"<<update->ntimestep<<",\"accepted_steps\":"<<accepted<<",\"termination_reason\":\""<<reasons[reason]<<"\",\"query_count\":"<<queries<<",\"invalid_query_count\":"<<invalid<<",\"nonfinite_count\":"<<nonfinite<<",\"kokkosable\":"<<kokkosable<<",\"datamask_read\":"<<datamask_read<<",\"datamask_modify\":"<<datamask_modify<<"}\n";}

#include "fix_frozen_flow_drag.hpp"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "update.h"
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <mpi.h>
using namespace LAMMPS_NS;
using namespace frozen;
FixFrozenFlowDrag::FixFrozenFlowDrag(LAMMPS*lmp,int narg,char**arg):Fix(lmp,narg,arg){
 if(narg!=12)error->all(FLERR,"frozen/flow/drag expects field mu stride prefix wall_at_spawn spawn_xyz margin");
 try{field=std::make_unique<FlowField>(arg[3]);mu=std::stod(arg[4]);stride=std::stoi(arg[5]);prefix=arg[6];wall_at_spawn=std::stod(arg[7]);for(int j=0;j<3;++j)spawn[j]=std::stod(arg[8+j]);margin=std::stod(arg[11]);if(mu<=0||stride<1)throw std::runtime_error("Invalid fix parameters");}catch(std::exception const&e){error->all(FLERR,e.what());}
 nevery=1;MPI_Comm_rank(world,&rank);trace.open(prefix+"_rank"+std::to_string(rank)+".csv");if(!trace)error->all(FLERR,"Cannot open coupling trajectory output");trace<<std::setprecision(17)<<"step,time_s,id,diameter_m,x,y,z,vx,vy,vz,ufx,ufy,ufz,Fx,Fy,Fz,applied_Fx,applied_Fy,applied_Fz,force_vx,force_vy,force_vz,center_wall_lower_bound_m,query_status\n";
 // Host fix, with standard conservative Fix data masks. ModifyKokkos supplies sync.
 // No custom CUDA kernels, no replacement of LAMMPS integration or drag equation.
}
int FixFrozenFlowDrag::setmask(){return FixConst::POST_FORCE|FixConst::END_OF_STEP|FixConst::POST_RUN;}
void FixFrozenFlowDrag::setup(int flag){post_force(flag);record();}
void FixFrozenFlowDrag::post_force(int){++calls;applied.clear();int bad=0;for(int i=0;i<atom->nlocal;++i)if(atom->mask[i]&groupbit){++queries;Vec x,v;for(int j=0;j<3;++j){x[j]=atom->x[i][j];v[j]=atom->v[i][j];}auto q=FlowFieldSampler(*field).query(x);if(q.status!=Status::VALID){++invalid;bad=1;continue;}try{Vec f=TechnicalStokesDrag::force(mu,2*atom->radius[i],q.u,v);for(int j=0;j<3;++j)atom->f[i][j]+=f[j];applied.emplace(atom->tag[i],Applied{f,v});}catch(...){++nonfinite;bad=1;}if(wall_at_spawn>=0){double ds=0;for(int j=0;j<3;++j)ds+=(x[j]-spawn[j])*(x[j]-spawn[j]);double lower=wall_at_spawn-std::sqrt(ds);minimum_clearance=std::min(minimum_clearance,lower);if(!(lower>atom->radius[i]+margin))bad=1;}}
 int global=0;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,world);if(global)error->all(FLERR,"INVALID_QUERY or nonfinite force or safe-region violation; no fallback");}
void FixFrozenFlowDrag::record(){int bad=0;bool write=update->ntimestep!=last_written&&(update->ntimestep%stride==0||update->ntimestep==update->laststep);for(int i=0;i<atom->nlocal;++i)if(atom->mask[i]&groupbit){Vec x,v;for(int j=0;j<3;++j){x[j]=atom->x[i][j];v[j]=atom->v[i][j];}auto q=FlowFieldSampler(*field).query(x);if(q.status!=Status::VALID){++invalid;bad=1;continue;}Vec f;try{f=TechnicalStokesDrag::force(mu,2*atom->radius[i],q.u,v);}catch(...){++nonfinite;bad=1;continue;}double lower=-1;if(wall_at_spawn>=0){double ds=0;for(int j=0;j<3;++j)ds+=(x[j]-spawn[j])*(x[j]-spawn[j]);lower=wall_at_spawn-std::sqrt(ds);minimum_clearance=std::min(minimum_clearance,lower);if(!(lower>atom->radius[i]+margin))bad=1;}if(write){double unwrapped[3];domain->unmap(atom->x[i],atom->image[i],unwrapped);auto it=applied.find(atom->tag[i]);if(it==applied.end()){bad=1;continue;}trace<<update->ntimestep<<','<<update->ntimestep*update->dt<<','<<atom->tag[i]<<','<<2*atom->radius[i];for(double a:unwrapped)trace<<','<<a;for(auto arr:{v,q.u,f,it->second.force,it->second.velocity})for(double a:arr)trace<<','<<a;trace<<','<<lower<<",0\n";}}
 int global=0;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,world);if(global)error->all(FLERR,"Completed-state coupling safety violation");if(write){last_written=update->ntimestep;if(update->ntimestep%100==0||update->ntimestep==update->laststep)trace.flush();}}
void FixFrozenFlowDrag::end_of_step(){record();}
void FixFrozenFlowDrag::post_run(){trace.flush();std::ofstream f(prefix+"_counters_rank"+std::to_string(rank)+".json");f<<std::setprecision(17)<<"{\"completed_step\":"<<update->ntimestep<<",\"post_force_calls\":"<<calls<<",\"force_queries\":"<<queries<<",\"invalid_query_count\":"<<invalid<<",\"nonfinite_count\":"<<nonfinite<<",\"host_fix\":true,\"kokkosable\":"<<kokkosable<<",\"datamask_read\":"<<datamask_read<<",\"datamask_modify\":"<<datamask_modify<<",\"min_center_wall_lower_bound_m\":"<<(wall_at_spawn>=0?minimum_clearance:-1)<<"}\n";}

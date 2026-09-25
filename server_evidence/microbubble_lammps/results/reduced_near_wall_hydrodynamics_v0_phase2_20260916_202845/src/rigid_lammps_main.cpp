#include "lammps.h"
#include "input.h"
#include "modify.h"
#include "fix.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "domain.h"
#include "update.h"
#include "library.h"
#include "rigid_math.hpp"
#include "flat_wall_engine.hpp"
#include <mpi.h>
#include <cmath>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <algorithm>
using namespace LAMMPS_NS;using namespace rigid;
static std::map<long long,std::array<double,9>>proposal;
class RigidStorage:public Fix{public:RigidStorage(LAMMPS*l,int n,char**a):Fix(l,n,a){time_integrate=1;nevery=1;}int setmask()override{return FixConst::INITIAL_INTEGRATE|FixConst::END_OF_STEP;}void initial_integrate(int)override{for(int i=0;i<atom->nlocal;i++){auto const&t=proposal.at(atom->tag[i]);for(int k=0;k<3;k++){atom->x[i][k]=t[k];atom->v[i][k]=t[k+3];atom->omega[i][k]=t[k+6];}}}void end_of_step()override{}};
static Fix*factory(LAMMPS*l,int n,char**a){return new RigidStorage(l,n,a);}
static std::vector<double>gather(const std::vector<double>&local,MPI_Comm world){int nr;MPI_Comm_size(world,&nr);int n=local.size();std::vector<int>c(nr),d(nr);MPI_Allgather(&n,1,MPI_INT,c.data(),1,MPI_INT,world);for(int k=1;k<nr;k++)d[k]=d[k-1]+c[k-1];std::vector<double>out(d.back()+c.back());MPI_Allgatherv(local.data(),n,MPI_DOUBLE,out.data(),c.data(),d.data(),MPI_DOUBLE,world);return out;}
static std::vector<Particle>particles(LAMMPS*l,int rank){std::vector<double>v;auto*a=l->atom;for(int i=0;i<a->nlocal;i++){v.push_back(a->tag[i]);v.push_back(rank);for(int k=0;k<3;k++)v.push_back(a->x[i][k]);v.push_back(a->radius[i]);}auto all=gather(v,l->world);std::vector<Particle>p;for(size_t k=0;k<all.size();k+=6)p.push_back({(long long)all[k],int(all[k+1]),{all[k+2],all[k+3],all[k+4]},all[k+5]});std::sort(p.begin(),p.end(),[](auto&a,auto&b){return a.id<b.id;});return p;}
static std::vector<std::pair<int,int>>neighbors(LAMMPS*l,const std::vector<Particle>&p,long long&crossrank){std::map<long long,int>idx;for(size_t i=0;i<p.size();i++)idx[p[i].id]=i;auto*list=l->force->pair->list;if(!list)throw std::runtime_error("LAMMPS_NEIGHBOR_LIST_UNAVAILABLE");std::vector<double>send;for(int k=0;k<list->inum;k++){int i=list->ilist[k];for(int jn=0;jn<list->numneigh[i];jn++){int j=list->firstneigh[i][jn]&NEIGHMASK;if(l->atom->tag[i]<l->atom->tag[j]){send.push_back(l->atom->tag[i]);send.push_back(l->atom->tag[j]);}}}auto v=gather(send,l->world);std::set<std::pair<int,int>>unique;for(size_t k=0;k<v.size();k+=2)unique.emplace(idx.at(v[k]),idx.at(v[k+1]));std::vector<std::pair<int,int>>pairs(unique.begin(),unique.end());for(auto ij:pairs){auto&a=p[ij.first];auto&b=p[ij.second];if(a.owner!=b.owner&&norm(sub(a.x,b.x))-a.a-b.a<.2*a.a*b.a/(a.a+b.a))crossrank++;}return pairs;}
static double swept(const std::vector<Particle>&p,const std::vector<Particle>&next,const std::vector<std::pair<int,int>>&pairs){double min=1e100;for(auto ij:pairs){int i=ij.first,j=ij.second;Vec r=sub(p[i].x,p[j].x),v=sub(sub(next[i].x,p[i].x),sub(next[j].x,p[j].x));double vv=dot(v,v),t=vv>0?std::clamp(-dot(r,v)/vv,0.,1.):0.;min=std::min(min,norm(add(r,mul(v,t)))-p[i].a-p[j].a);}return min;}
int main(int argc,char**argv) {
 MPI_Init(&argc,&argv);int rank,nr;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nr);
 LAMMPS*l=nullptr;
 try {
  if(argc!=2)throw std::runtime_error("Usage: rigid_lmp_cpu CASE_CONFIG; CPU analytic plane only");
  std::ifstream in(argv[1]);if(!in)throw std::runtime_error("CONFIG_OPEN_FAILED");
  std::map<std::string,std::string>cfg;std::string key,value;while(in>>key>>value)cfg[key]=value;
  auto get=[&](std::string k){return cfg.at(k);};
  auto number=[&](std::string k){double v=std::stod(get(k));if(!std::isfinite(v))throw std::runtime_error("NONFINITE_CONFIG");return v;};
  auto vec=[&](std::string p){return Vec{number(p+"x"),number(p+"y"),number(p+"z")};};
  if(get("shear_source")!="ANALYTIC"||get("geometry")!="ANALYTIC_PLANE")throw std::runtime_error("UNSUPPORTED_SOURCE");
  const double dtmax=number("dt_max"),tmax=number("max_time");const int stride=std::stoi(get("stride")),stepscap=std::stoi(get("steps_cap"));
  if(dtmax<=0||tmax<0||stride<1||stepscap<0)throw std::runtime_error("INVALID_TIME_CONFIG");
  const bool dumpmatrix=get("dump_matrix")=="YES";
  phase2::Plane plane{vec("origin_"),{vec("t1_"),vec("t2_"),vec("n_")},number("gamma_dot")};
  phase2::FlatWallEngine wall(get("resistance_table"),plane);
  std::vector<std::string>args={"rigid_lmp_cpu","-screen","none","-log","none"};std::vector<char*>av;for(auto&s:args)av.push_back(s.data());
  l=new LAMMPS(av.size(),av.data(),MPI_COMM_WORLD);(*l->modify->fix_map)["sonovue/rigid/storage"]=&factory;
  auto command=[&](std::string s){l->input->one(s);};
  command("units si");command("atom_style sphere");command("boundary f f f");command("newton off");command("atom_modify map array");
  command("processors "+std::to_string(nr)+" 1 1");command("read_data "+get("data"));
  command("pair_style zero 8e-6 full");command("pair_coeff * *");command("neighbor 4e-7 bin");command("neigh_modify every 1 delay 0 check no");command("comm_modify vel yes");
  command("fix rigid all sonovue/rigid/storage");command("thermo 1000000");command("run 0 post no");
  auto p=particles(l,rank);if(p.size()!=1)throw std::runtime_error("PHASE2_REQUIRES_SINGLE_PARTICLE");
  long long crossrank=0;auto pairs=neighbors(l,p,crossrank);if(!pairs.empty())throw std::runtime_error("PHASE2_PAIR_LIST_NOT_EMPTY");
  double time=0,dt=0,maxres=0,min_swept=1e100,max_storage_error=0;int accepted=0,ownerchanges=0;long long evaluations=0;
  std::ofstream trajectory,stages,matrices;
  if(rank==0) {
   trajectory.open("TRAJECTORIES.csv");stages.open("INTEGRATION_STAGES.csv");
   trajectory<<std::setprecision(17);stages<<std::setprecision(17);
   trajectory<<"time,step,particle_id,x,y,z,radius,gap,epsilon,Vx,Vy,Vz,Omega_x,Omega_y,Omega_z,FU,FOMEGA,normal_velocity,tangential_velocity,wall_status,kernel_status,owner_rank,linear_residual\n";
   stages<<"evaluation_id,step,stage,time,dt,particle_id,x,y,z,gap,epsilon,FU,FOMEGA,Vx,Vy,Vz,Omega_x,Omega_y,Omega_z,normal_velocity,wall_status,kernel_status,linear_residual,pair_count\n";
   if(dumpmatrix) {
    matrices.open("MATRIX_RHS.csv");matrices<<std::setprecision(17)<<"evaluation_id,stage,epsilon";
    for(int i=0;i<36;++i)matrices<<",Rwall_"<<i;
    for(int i=0;i<36;++i)matrices<<",Rtotal_"<<i;
    for(int i=0;i<6;++i)matrices<<",bwall_"<<i;
    for(int i=0;i<6;++i)matrices<<",rhs_"<<i;
    for(int i=0;i<6;++i)matrices<<",q_"<<i;
    matrices<<'\n';
   }
  }
  auto background=[&](const std::vector<Particle>&pp){std::vector<Background>b;for(auto const&x:pp)b.push_back(wall.background(x));return b;};
  std::vector<reducedwall::ReducedWallShearContribution>last_walls;
  // Every invocation queries the supplied stage positions, including the midpoint.
  auto hydro=[&](const std::vector<Particle>&pp,const std::vector<Background>&bg,
                const std::vector<Particle>*base,double subdt,double t,const char*stage) {
   std::vector<reducedwall::ReducedWallShearContribution>w;
   for(auto const&x:pp){auto c=wall.evaluate(x);if(c.status!=reducedwall::Status::OK)throw std::runtime_error(phase2::engine_status(c.status));w.push_back(c);}
   auto sol=solve(pp,bg,pairs,base,subdt,&w);++evaluations;maxres=std::max(maxres,sol.residual);
   if(rank==0)for(size_t i=0;i<pp.size();++i) {
    stages<<evaluations<<','<<accepted<<','<<stage<<','<<t<<','<<subdt<<','<<pp[i].id;
    for(double x:pp[i].x)stages<<','<<x;
    stages<<','<<phase2::gap(plane,pp[i].x,pp[i].a)<<','<<w[i].resistance.epsilon<<','<<w[i].FU<<','<<w[i].FOMEGA;
    for(int k=0;k<6;++k)stages<<','<<sol.q[6*i+k];
    stages<<','<<dot({sol.q[6*i],sol.q[6*i+1],sol.q[6*i+2]},plane.frame.n)<<",OK,OK,"<<sol.residual<<','<<pairs.size()<<'\n';
    if(dumpmatrix) {
     auto sys=assemble(pp,bg,pairs,&w);matrices<<evaluations<<','<<stage<<','<<w[i].resistance.epsilon;
     for(double v:w[i].resistance.wall_excess_SI)matrices<<','<<v;
     for(double v:sys.R)matrices<<','<<v;
     for(double v:w[i].b_wall_shear)matrices<<','<<v;
     for(double v:sys.rhs)matrices<<','<<v;
     for(double v:sol.q)matrices<<','<<v;
     matrices<<'\n';
    }
   }
   last_walls=std::move(w);return sol;
  };
  auto record=[&](const Solution&s) {
   if(rank)return;
   for(size_t i=0;i<p.size();++i) {
    const auto&w=last_walls[i];trajectory<<time<<','<<accepted<<','<<p[i].id;
    for(double x:p[i].x)trajectory<<','<<x;
    trajectory<<','<<p[i].a<<','<<phase2::gap(plane,p[i].x,p[i].a)<<','<<w.resistance.epsilon;
    for(int k=0;k<6;++k)trajectory<<','<<s.q[6*i+k];
    Vec v{s.q[6*i],s.q[6*i+1],s.q[6*i+2]};
    trajectory<<','<<w.FU<<','<<w.FOMEGA<<','<<dot(v,plane.frame.n)<<','<<dot(v,plane.frame.t1)<<",OK,OK,"<<p[i].owner<<','<<s.residual<<'\n';
   }
  };
  auto bg=background(p);auto current=hydro(p,bg,nullptr,0,time,"INITIAL");record(current);
  while(time<tmax-1e-16&&accepted<stepscap) {
   dt=std::min(dtmax,tmax-time);
   // Inherited production displacement limit for the padded LAMMPS neighbor list.
   for(size_t i=0;i<p.size();++i){double v=norm({current.q[6*i],current.q[6*i+1],current.q[6*i+2]});if(v>0)dt=std::min(dt,.25*1.9989918081065344e-7/v);}
   auto first=hydro(p,bg,&p,.5*dt,time,"START");auto mid=p;
   for(size_t i=0;i<p.size();++i)for(int k=0;k<3;++k)mid[i].x[k]+=.5*dt*first.q[6*i+k];
   for(size_t i=0;i<p.size();++i)if(!phase2::segment_safe(plane,p[i].x,mid[i].x,p[i].a))throw std::runtime_error("HARD_WALL_PROPOSAL_REJECTED");
   auto mbg=background(mid);auto second=hydro(mid,mbg,&p,dt,time+.5*dt,"MIDPOINT");auto next=p;
   for(size_t i=0;i<p.size();++i)for(int k=0;k<3;++k)next[i].x[k]+=dt*second.q[6*i+k];
   for(size_t i=0;i<p.size();++i)if(!phase2::segment_safe(plane,p[i].x,next[i].x,p[i].a))throw std::runtime_error("HARD_WALL_PROPOSAL_REJECTED");
   const double sg=swept(p,next,pairs);if(sg< -1e-12)throw std::runtime_error("SWEPT_NONOVERLAP_GATE_FAILED");min_swept=std::min(min_swept,sg);
   auto nbg=background(next);auto finalq=hydro(next,nbg,&next,dt,time+dt,"ENDPOINT");
   proposal.clear();for(size_t i=0;i<p.size();++i){std::array<double,9>v;for(int k=0;k<3;++k){v[k]=next[i].x[k];v[k+3]=finalq.q[6*i+k];v[k+6]=finalq.q[6*i+k+3];}proposal[p[i].id]=v;}
   std::ostringstream ts;ts<<std::setprecision(17)<<"timestep "<<dt;command(ts.str());command("run 1 pre no post no");
   double storage_error=0;for(int i=0;i<l->atom->nlocal;++i){const auto&v=proposal.at(l->atom->tag[i]);for(int k=0;k<3;++k){storage_error=std::max(storage_error,std::abs(l->atom->x[i][k]-v[k]));storage_error=std::max(storage_error,std::abs(l->atom->v[i][k]-v[k+3]));storage_error=std::max(storage_error,std::abs(l->atom->omega[i][k]-v[k+6]));}}
   double allerror=0;MPI_Allreduce(&storage_error,&allerror,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);max_storage_error=std::max(max_storage_error,allerror);if(allerror!=0)throw std::runtime_error("LAMMPS_STORAGE_MISMATCH");
   const int oldowner=p.front().owner;accepted++;time+=dt;p=particles(l,rank);if(p.size()!=1)throw std::runtime_error("PARTICLE_LOST_OR_DUPLICATED");ownerchanges+=p.front().owner!=oldowner;
   pairs=neighbors(l,p,crossrank);if(!pairs.empty())throw std::runtime_error("PHASE2_PAIR_LIST_NOT_EMPTY");
   bg=background(p);current=hydro(p,bg,&p,dt,time,"ACCEPTED");
   if(accepted%stride==0||time>=tmax-1e-16)record(current);
  }
  if(time<tmax-1e-16)throw std::runtime_error("STEP_CAP_NOT_COMPLETE");
  if(rank==0){std::ofstream s("RUN_STATE.json");s<<std::setprecision(17)<<"{\"status\":\"COMPLETE\",\"reason\":\"MAX_PHYSICAL_TIME\",\"accepted_steps\":"<<accepted<<",\"time_s\":"<<time<<",\"max_residual\":"<<maxres<<",\"mpi_ranks\":"<<nr<<",\"ownership_changes\":"<<ownerchanges<<",\"hydrodynamic_evaluations\":"<<evaluations<<",\"max_LAMMPS_storage_error\":"<<max_storage_error<<",\"pair_count\":0,\"flow_source\":\"ANALYTIC\",\"physical_integrator\":\"stable midpoint RK2 with LAMMPS storage/migration\",\"wall_domain\":[0.001,0.2],\"flow_feedback\":false}\n";
   std::cout<<"COMPLETE steps="<<accepted<<" evaluations="<<evaluations<<" owner_changes="<<ownerchanges<<std::endl;}
  delete l;l=nullptr;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();
 }catch(std::exception const&e){if(rank==0){std::ofstream s("RUN_STATE.json");s<<"{\"status\":\"REJECTED\",\"reason\":\""<<e.what()<<"\"}\n";std::cerr<<"RIGID_ENGINE_REJECTED "<<e.what()<<std::endl;}MPI_Abort(MPI_COMM_WORLD,2);}
 MPI_Finalize();return 0;
}

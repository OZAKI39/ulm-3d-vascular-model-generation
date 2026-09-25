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
#include "wall_distance.hpp"
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
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rank,nr;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nr);int code=0;LAMMPS*l=nullptr;try{
 if(argc<2)throw std::runtime_error("Usage rigid_lmp CASE_CONFIG [kokkos]");std::ifstream in(argv[1]);std::map<std::string,std::string>cfg;std::string key,value;while(in>>key>>value)cfg[key]=value;
 auto get=[&](std::string k){return cfg.at(k);};double dtmax=std::stod(get("dt_max")),c_gap=std::stod(get("c_gap")),tmax=std::stod(get("max_time")),margin=std::stod(get("margin"));int stride=std::stoi(get("stride")),mode=std::stoi(get("drive_mode")),stepscap=std::stoi(get("steps_cap"));double drivescale=std::stod(get("drive_scale"));
 frozen::FlowField field(get("field"));FlowGradientSampler sampler(field);std::unique_ptr<passive::WallDistance>wall;if(get("wall")!="NONE")wall=std::make_unique<passive::WallDistance>(get("wall"));
 std::vector<std::string>args={"rigid_lmp","-screen","none","-log","none"};bool gpu=argc>2;if(gpu){for(auto s:{"-k","on","g","1","-sf","kk","-pk","kokkos","neigh","half","newton","off","gpu/aware","off"})args.emplace_back(s);}std::vector<char*>av;for(auto&s:args)av.push_back(s.data());l=new LAMMPS(av.size(),av.data(),MPI_COMM_WORLD);(*l->modify->fix_map)["sonovue/rigid/storage"]=&factory;
 auto command=[&](std::string s){l->input->one(s);};command("units si");command("atom_style sphere");command("boundary f f f");command("newton off");command("atom_modify map array");command("read_data "+get("data"));command("pair_style zero 8e-6 full");command("pair_coeff * *");command("neighbor 4e-7 bin");command("neigh_modify every 1 delay 0 check no");command("comm_modify vel yes");command("fix rigid all sonovue/rigid/storage");command("thermo 1000000");command("run 0 post no");
 auto p=particles(l,rank);long long crossrank=0;auto pairs=neighbors(l,p,crossrank);double time=0,dt=0,min_swept=1e100,maxres=0;int accepted=0,overlaps=0,constraint_steps=0,gap_limited=0;std::string reason="MAX_PHYSICAL_TIME";
 auto background=[&](const std::vector<Particle>&p){std::vector<Background>b;for(auto const&x:p){auto z=sampler.query(x.x);if(mode==1)z.dw[2]+=drivescale;if(mode==2)z.dw[2]+=x.id%2?drivescale:-drivescale;if(mode==3){for(int k=0;k<3;k++){z.du[k]+=drivescale*std::sin(x.id*1.7+k);z.dw[k]+=1e6*drivescale*std::cos(x.id*2.3+k);}}b.push_back(z);}return b;};
 std::ofstream trajectory,events,history,stages;if(rank==0){trajectory.open("TRAJECTORIES.csv");events.open("PAIR_HYDRODYNAMIC_EVENTS.csv");history.open("SOLVER_HISTORY.csv");stages.open("INTEGRATION_STAGES.csv");for(auto*f:{&trajectory,&events,&history,&stages})*f<<std::setprecision(17);
 trajectory<<"particle_id,time_s,step,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,omega_x_rad_s,omega_y_rad_s,omega_z_rad_s,diameter_um,angular_speed_rad_s,nearest_gap_m,active_pair_count,center_wall_distance_m,surface_wall_gap_m,owner_rank\n";
 events<<"time_s,step,id_i,id_j,diameter_i_um,diameter_j_um,gap_m,normal_relative_speed_m_s,tangential_relative_speed_m_s,omega_ix,omega_iy,omega_iz,omega_jx,omega_jy,omega_jz,normal_resistance_kg_s,shear_resistance_kg_s,pump_resistance_kg_m2_s,twist_resistance,force_ix,force_iy,force_iz,force_jx,force_jy,force_jz,torque_ix,torque_iy,torque_iz,torque_jx,torque_jy,torque_jz,normal_force_x,normal_force_y,normal_force_z,shear_force_x,shear_force_y,shear_force_z,pump_torque_ix,pump_torque_iy,pump_torque_iz,twist_contribution,constraint_active\n";
 history<<"step,time_s,dt_s,max_residual,min_swept_gap_m,constraints,cross_rank_active_pair_observations,neighbor_pairs\n";stages<<"step,stage,time_s,dt_s,particle_id,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,omega_x,omega_y,omega_z\n";}
 auto stage_record=[&](int stage,const std::vector<Particle>&pp,const Solution&s,double t,double dt){if(rank)return;for(size_t i=0;i<pp.size();i++){stages<<accepted<<','<<stage<<','<<t<<','<<dt<<','<<pp[i].id;for(double x:pp[i].x)stages<<','<<x;for(int k=0;k<6;k++)stages<<','<<s.q[6*i+k];stages<<'\n';}};
 int last_recorded=-1;auto record=[&](const Solution&s,bool force){if(rank||last_recorded==accepted)return;if(!force&&accepted%stride)return;for(size_t i=0;i<p.size();i++){double gap=1e100;int active=0;for(auto ij:pairs)if(ij.first==int(i)||ij.second==int(i)){int j=ij.first==int(i)?ij.second:ij.first;double h=norm(sub(p[i].x,p[j].x))-p[i].a-p[j].a;gap=std::min(gap,h);if(h<.2*p[i].a*p[j].a/(p[i].a+p[j].a))active++;}double wd=wall?wall->distance(p[i].x):-1;trajectory<<p[i].id<<','<<time<<','<<accepted;for(double x:p[i].x)trajectory<<','<<x;for(int k=0;k<6;k++)trajectory<<','<<s.q[6*i+k];trajectory<<','<<2e6*p[i].a<<','<<norm({s.q[6*i+3],s.q[6*i+4],s.q[6*i+5]})<<','<<gap<<','<<active<<','<<wd<<','<<(wall?wd-p[i].a:-1.)<<','<<p[i].owner<<'\n';}
  for(auto const&e:s.events){Vec v{};for(int k=0;k<3;k++)v[k]=s.q[6*e.i+k]-s.q[6*e.j+k];double vn=dot(v,e.n);events<<time<<','<<accepted<<','<<p[e.i].id<<','<<p[e.j].id<<','<<2e6*p[e.i].a<<','<<2e6*p[e.j].a<<','<<e.gap<<','<<vn<<','<<norm(sub(v,mul(e.n,vn)));for(int i:{e.i,e.j})for(int k=3;k<6;k++)events<<','<<s.q[6*i+k];events<<','<<e.c.sq<<','<<e.c.sh<<','<<e.c.pu<<",PENDING_NOT_IMPLEMENTED";for(Vec v:{e.force,mul(e.force,-1),e.ti,e.tj,e.fn,e.fs,e.tpu})for(double x:v)events<<','<<x;events<<",NOT_IMPLEMENTED,"<<e.constraint<<'\n';}trajectory.flush();events.flush();last_recorded=accepted;};
 auto bg=background(p);auto current=solve(p,bg,pairs);record(current,true);
 while(time<tmax-1e-16&&accepted<stepscap){
  dt=std::min(dtmax,tmax-time);for(auto const&e:current.events){Vec v;for(int k=0;k<3;k++)v[k]=current.q[6*e.i+k]-current.q[6*e.j+k];double vn=std::abs(dot(v,e.n));double re=p[e.i].a*p[e.j].a/(p[e.i].a+p[e.j].a);double candidate=c_gap*std::max(e.gap,.001*re)/std::max(vn,1e-30);if(candidate<dt){dt=candidate;gap_limited++;}}
  // Reject a too-large displacement before it can outrun the padded LAMMPS list.
  for(size_t i=0;i<p.size();i++){double v=norm({current.q[6*i],current.q[6*i+1],current.q[6*i+2]});if(v>0)dt=std::min(dt,.25*1.9989918081065344e-7/v);}
  auto first=solve(p,bg,pairs,&p,.5*dt);auto mid=p;for(size_t i=0;i<p.size();i++)for(int k=0;k<3;k++)mid[i].x[k]+=.5*dt*first.q[6*i+k];bool safe=true;for(size_t i=0;i<p.size();i++)if(wall&&!wall->segment_safe(p[i].x,mid[i].x,p[i].a+margin))safe=false;if(!safe){reason="WALL_SAFETY_STOP";break;}
  auto mbg=background(mid);auto second=solve(mid,mbg,pairs,&p,dt);auto next=p;for(size_t i=0;i<p.size();i++)for(int k=0;k<3;k++)next[i].x[k]+=dt*second.q[6*i+k];for(size_t i=0;i<p.size();i++)if(wall&&!wall->segment_safe(p[i].x,next[i].x,p[i].a+margin))safe=false;if(!safe){reason="WALL_SAFETY_STOP";break;}
  double sg=swept(p,next,pairs);if(sg < -1e-12)throw std::runtime_error("SWEPT_NONOVERLAP_GATE_FAILED");min_swept=std::min(min_swept,sg);if(sg < -1e-12)overlaps++;auto nbg=background(next);auto finalq=solve(next,nbg,pairs,&next,dt);maxres=std::max({maxres,first.residual,second.residual,finalq.residual});constraint_steps+=!second.constraints.empty();stage_record(0,p,first,time,dt);stage_record(1,mid,second,time+.5*dt,dt);
  proposal.clear();for(size_t i=0;i<p.size();i++){std::array<double,9>v;for(int k=0;k<3;k++){v[k]=next[i].x[k];v[k+3]=finalq.q[6*i+k];v[k+6]=finalq.q[6*i+k+3];}proposal[p[i].id]=v;}
  std::ostringstream ts;ts<<std::setprecision(17)<<"timestep "<<dt;command(ts.str());command("run 1 pre no post no");accepted++;time+=dt;p=particles(l,rank);pairs=neighbors(l,p,crossrank);bg=background(p);current=solve(p,bg,pairs,&p,dt);record(current,time>=tmax-1e-16);if(rank==0){history<<accepted<<','<<time<<','<<dt<<','<<maxres<<','<<sg<<','<<second.constraints.size()<<','<<crossrank<<','<<pairs.size()<<'\n';if(accepted%100==0){history.flush();stages.flush();std::cout<<"accepted="<<accepted<<" time="<<time<<" pairs="<<current.events.size()<<std::endl;}}
 }
 if(accepted>=stepscap&&time<tmax-1e-16)reason="STEP_CAP_NOT_COMPLETE";if(reason!="MAX_PHYSICAL_TIME")record(current,true);
 if(rank==0){std::ofstream s("RUN_STATE.json");s<<std::setprecision(17)<<"{\"status\":\"TERMINAL\",\"reason\":\""<<reason<<"\",\"accepted_steps\":"<<accepted<<",\"time_s\":"<<time<<",\"max_residual\":"<<maxres<<",\"min_swept_gap_m\":"<<min_swept<<",\"accepted_overlap_count\":"<<overlaps<<",\"constraint_steps\":"<<constraint_steps<<",\"gap_limited_decisions\":"<<gap_limited<<",\"cross_rank_active_pair_observations\":"<<crossrank<<",\"mpi_ranks\":"<<nr<<",\"kokkos_compatibility\":"<<(gpu?"true":"false")<<",\"neighbor_source\":\"LAMMPS pair_zero full list\",\"O_N2_particle_search\":false,\"physical_integrator\":\"project RK2 with LAMMPS storage/migration\"}\n";std::cout<<"TERMINAL "<<reason<<" time="<<time<<std::endl;}
 delete l;l=nullptr;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();
 }catch(std::exception const&e){std::cerr<<"RIGID_ENGINE_FAILED "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);code=2;}MPI_Finalize();return code;}

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
#include "wall_model.hpp"
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
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,nr;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nr);LAMMPS*l=nullptr;
 try{
 if(argc<2)throw std::runtime_error("Usage wall_lmp case.cfg [kokkos]");std::ifstream input(argv[1]);std::map<std::string,std::string>cfg;std::string key,value;while(input>>key>>value)cfg[key]=value;
 auto get=[&](std::string k){return cfg.at(k);};double dtmax=std::stod(get("dt_max")),cgap=std::stod(get("c_gap")),cwall=std::stod(get("c_wall")),tmax=std::stod(get("max_time"));int stride=std::stoi(get("stride")),cap=std::stoi(get("steps_cap"));
 bool hydro=get("wall_resistance")=="1",hard=get("hard_wall")=="1",pair_on=get("pair_enabled")=="1",control=get("control_overlap")=="1";int pattern=std::stoi(get("drive_pattern"));Vec drive={std::stod(get("drive_x")),std::stod(get("drive_y")),std::stod(get("drive_z"))},torque={std::stod(get("torque_x")),std::stod(get("torque_y")),std::stod(get("torque_z"))};
 frozen::FlowField field(get("field"));FlowGradientSampler sampler(field);wallv0::WallResistanceLookup lookup(get("table"));wallv0::STLNearestWallQuery geometry(get("wall"));
 std::vector<std::string>args={"wall_lmp","-screen","none","-log","none"};bool gpu=argc>2;if(gpu)for(auto s:{"-k","on","g","1","-sf","kk","-pk","kokkos","neigh","half","newton","off","gpu/aware","off"})args.emplace_back(s);std::vector<char*>av;for(auto&s:args)av.push_back(s.data());l=new LAMMPS(av.size(),av.data(),MPI_COMM_WORLD);(*l->modify->fix_map)["sonovue/wall/storage"]=&factory;
 auto command=[&](std::string s){l->input->one(s);};command("units si");command("atom_style sphere");command("boundary f f f");command("newton off");command("atom_modify map array");command("read_data "+get("data"));command("pair_style zero 8e-6 full");command("pair_coeff * *");command("neighbor 4e-7 bin");command("neigh_modify every 1 delay 0 check no");command("comm_modify vel yes");command("fix wall all sonovue/wall/storage");command("thermo 1000000");command("run 0 post no");
 auto p=particles(l,rank);int expected=p.size();long long crossrank=0;auto allpairs=neighbors(l,p,crossrank);auto pairs=pair_on?allpairs:std::vector<std::pair<int,int>>{};
 std::map<long long,wallv0::Frame>previous;double time=0,min_wall=1e100,min_pair=1e100,maxres=0;long long ownership_changes=0,wall_activations=0,wall_active_observations=0;int accepted=0,rejected=0,invalid=0,wall_overlap=0,pair_overlap=0,nonfinite=0,wall_limited=0;std::string reason="MAX_PHYSICAL_TIME";
 auto background=[&](const std::vector<Particle>&pp){std::vector<Background>b;for(auto const&x:pp){auto z=sampler.query(x.x);if(pattern==0||x.id%3==0){z.du=drive;z.dw=torque;}b.push_back(z);}return b;};
 auto inspect=[&](const std::vector<Particle>&pp){std::vector<wallv0::Entry>out;for(auto const&x:pp){wallv0::Entry w;w.geometry=geometry.query(x.x,x.a);if(control&&get("wall")=="FLAT")w.geometry.valid=w.geometry.planar=true;auto prev=previous.find(x.id);w.Q=wallv0::frame(w.geometry.normal,prev==previous.end()?nullptr:&prev->second);w.validity=lookup.validity_class(w.geometry.gap/x.a);w.active=hydro&&w.geometry.gap/x.a<=20;w.excess=(w.geometry.planar||!w.active)?lookup.global_excess(x.a,.001,w.geometry.gap,w.Q,hydro):wallv0::Matrix6{};out.push_back(w);}return out;};
 auto check=[&](const std::vector<wallv0::Entry>&entries){for(auto const&w:entries){if(!w.geometry.valid){invalid++;return std::string("INVALID_QUERY");}if(w.active&&!w.geometry.planar)return std::string("BLOCKED_LOCAL_PLANE_VALIDITY");}return std::string();};
 std::ofstream trajectory,events,history,stages;if(rank==0){trajectory.open("TRAJECTORIES.csv");events.open("WALL_INTERACTION_EVENTS.csv");history.open("SOLVER_HISTORY.csv");stages.open("INTEGRATION_STAGES.csv");for(auto*f:{&trajectory,&events,&history,&stages})*f<<std::setprecision(17);
 trajectory<<"particle_id,time_s,step,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,omega_x_rad_s,omega_y_rad_s,omega_z_rad_s,diameter_um,radius_m,nearest_x_m,nearest_y_m,nearest_z_m,center_wall_distance_m,surface_wall_gap_m,epsilon,normal_x,normal_y,normal_z,triangle_id,rms_over_a,normal_spread_deg,normal_angle_deg,table_validity_class,wall_active,constraint_active,owner_rank,local_plane_valid,wall_excess_normal_scaled,wall_excess_parallel_scaled,wall_excess_RR_parallel_scaled,wall_excess_RR_normal_scaled,wall_excess_TR_scaled\n";
 events<<"time_s,particle_id,gap_um,epsilon,normal_velocity_m_s,tangential_speed_m_s,angular_speed_rad_s,normal_wall_resistance_kg_s,parallel_wall_resistance_kg_s,rotational_parallel_kg_m2_s,rotational_normal_kg_m2_s,TR_kg_m_s,constraint_active,triangle_id,rms_over_a,normal_spread_deg\n";
 history<<"step,time_s,dt_s,max_residual,min_swept_wall_gap_bound_m,min_swept_pair_gap_m,wall_constraints,pair_constraints,rejected_trials,ownership_changes\n";
 stages<<"step,stage,time_s,dt_s,particle_id,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,omega_x_rad_s,omega_y_rad_s,omega_z_rad_s\n";}
 auto stage_record=[&](int stage,const std::vector<Particle>&pp,const Solution&s,double tt,double dt){if(rank)return;for(size_t i=0;i<pp.size();i++){stages<<accepted<<','<<stage<<','<<tt<<','<<dt<<','<<pp[i].id;for(double v:pp[i].x)stages<<','<<v;for(int j=0;j<6;j++)stages<<','<<s.q[6*i+j];stages<<'\n';}};
 int last=-1;auto record=[&](const Solution&s,const std::vector<wallv0::Entry>&ww,bool force){if(rank||last==accepted||(!force&&accepted%stride))return;for(size_t i=0;i<p.size();i++){
  auto const&g=ww[i].geometry;bool constrained=std::find(s.wall_constraints.begin(),s.wall_constraints.end(),i)!=s.wall_constraints.end();auto E=lookup.scaled_excess(g.gap/p[i].a);if(!hydro)E={};
  trajectory<<p[i].id<<','<<time<<','<<accepted;for(double v:p[i].x)trajectory<<','<<v;for(int j=0;j<6;j++)trajectory<<','<<s.q[6*i+j];trajectory<<','<<2e6*p[i].a<<','<<p[i].a;for(double v:g.closest)trajectory<<','<<v;trajectory<<','<<g.distance<<','<<g.gap<<','<<g.gap/p[i].a;for(double v:g.normal)trajectory<<','<<v;trajectory<<','<<g.triangle_id<<','<<g.rms_over_a<<','<<g.normal_spread_deg<<','<<g.nearest_normal_angle_deg<<','<<ww[i].validity<<','<<ww[i].active<<','<<constrained<<','<<p[i].owner<<','<<g.planar<<','<<E[14]<<','<<E[0]<<','<<E[21]<<','<<E[35]<<','<<E[4]<<'\n';
  if(ww[i].active||constrained){Vec v={s.q[6*i],s.q[6*i+1],s.q[6*i+2]},omega={s.q[6*i+3],s.q[6*i+4],s.q[6*i+5]};double vn=dot(g.normal,v),bt=6*std::acos(-1.)*.001*p[i].a,br=8*std::acos(-1.)*.001*std::pow(p[i].a,3);events<<time<<','<<p[i].id<<','<<g.gap*1e6<<','<<g.gap/p[i].a<<','<<vn<<','<<norm(sub(v,mul(g.normal,vn)))<<','<<norm(omega)<<','<<E[14]*bt<<','<<E[0]*bt<<','<<E[21]*br<<','<<E[35]*br<<','<<E[4]*std::sqrt(bt*br)<<','<<constrained<<','<<g.triangle_id<<','<<g.rms_over_a<<','<<g.normal_spread_deg<<'\n';}
 }trajectory.flush();events.flush();last=accepted;};
 auto ww=inspect(p);auto blocked=check(ww);auto bg=background(p);auto wb=wallv0::blocks(ww,ww,hard);Solution current;
 if(blocked.empty())current=solve(p,bg,pairs,nullptr,0,&wb);else{reason=blocked;current.q.assign(6*p.size(),0);}
 record(current,ww,true);for(auto const&w:ww)min_wall=std::min(min_wall,w.geometry.gap);
 while(reason=="MAX_PHYSICAL_TIME"&&time<tmax-1e-16&&accepted<cap){
  double dt=std::min(dtmax,tmax-time);
  for(auto const&e:current.events){Vec v;for(int j=0;j<3;j++)v[j]=current.q[6*e.i+j]-current.q[6*e.j+j];double speed=std::abs(dot(v,e.n)),re=p[e.i].a*p[e.j].a/(p[e.i].a+p[e.j].a);dt=std::min(dt,cgap*std::max(e.gap,.001*re)/std::max(speed,1e-30));}
  for(size_t i=0;i<p.size();i++){Vec v={current.q[6*i],current.q[6*i+1],current.q[6*i+2]};double speed=norm(v);if(speed>0)dt=std::min(dt,.25*field.dx/speed);double vn=dot(ww[i].geometry.normal,v);if((hydro||hard)&&vn<0){double limit=cwall*std::max(ww[i].geometry.gap,.001*p[i].a)/(-vn);if(limit<dt){dt=limit;wall_limited++;}}}
  bool accepted_trial=false;std::vector<Particle>mid,next;std::vector<wallv0::Entry>nw;Solution first,second,finalq;double sg=1e100;
  for(int trial=0;trial<40;trial++){
   if(dt<1e-15){reason="NO_SAFE_SUBSTEP";break;}
   first=solve(p,bg,pairs,&p,.5*dt,&wb);mid=p;for(size_t i=0;i<p.size();i++)for(int j=0;j<3;j++)mid[i].x[j]+=.5*dt*first.q[6*i+j];bool safe=true;
   if(hard)for(size_t i=0;i<p.size();i++)if(!geometry.swept_safe(p[i].x,mid[i].x,p[i].a))safe=false;
   if(swept(p,mid,allpairs)<-1e-12)safe=false;
   if(!safe){dt*=.5;rejected++;continue;}
   auto mw=inspect(mid);blocked=check(mw);if(!blocked.empty()){reason=blocked;break;}
   auto mb=wallv0::blocks(mw,ww,hard);auto mbg=background(mid);second=solve(mid,mbg,pairs,&p,dt,&mb);next=p;for(size_t i=0;i<p.size();i++)for(int j=0;j<3;j++)next[i].x[j]+=dt*second.q[6*i+j];
   if(hard)for(size_t i=0;i<p.size();i++)if(!geometry.swept_safe(p[i].x,next[i].x,p[i].a))safe=false;
   sg=swept(p,next,allpairs);if(sg<-1e-12)safe=false;
   if(!safe){dt*=.5;rejected++;continue;}
   nw=inspect(next);blocked=check(nw);if(!blocked.empty()){reason=blocked;break;}
   auto nb=wallv0::blocks(nw,nw,hard);auto nbg=background(next);finalq=solve(next,nbg,pairs,&next,dt,&nb);accepted_trial=true;break;
  }
  if(!accepted_trial){if(reason=="MAX_PHYSICAL_TIME")reason="NO_SAFE_SUBSTEP";break;}
  maxres=std::max({maxres,first.residual,second.residual,finalq.residual});wall_activations+=second.wall_constraints.size();min_pair=std::min(min_pair,sg);stage_record(0,p,first,time,dt);stage_record(1,mid,second,time+.5*dt,dt);
  for(size_t i=0;i<p.size();i++){min_wall=std::min(min_wall,nw[i].geometry.gap);if(nw[i].geometry.gap<-1e-12)wall_overlap++;wall_active_observations+=nw[i].active;previous[p[i].id]=nw[i].Q;}if(sg<-1e-12)pair_overlap++;
  proposal.clear();for(size_t i=0;i<p.size();i++){std::array<double,9>v;for(int j=0;j<3;j++){v[j]=next[i].x[j];v[j+3]=finalq.q[6*i+j];v[j+6]=finalq.q[6*i+j+3];if(!std::isfinite(v[j])||!std::isfinite(v[j+3])||!std::isfinite(v[j+6]))nonfinite++;}proposal[p[i].id]=v;}if(nonfinite)throw std::runtime_error("NONFINITE_STATE");
  std::ostringstream ts;ts<<std::setprecision(17)<<"timestep "<<dt;command(ts.str());command("run 1 pre no post no");accepted++;time+=dt;auto oldp=p;p=particles(l,rank);if(int(p.size())!=expected)throw std::runtime_error("LOST_ATOMS");for(size_t i=0;i<p.size();i++)ownership_changes+=p[i].owner!=oldp[i].owner;
  allpairs=neighbors(l,p,crossrank);pairs=pair_on?allpairs:std::vector<std::pair<int,int>>{};ww=inspect(p);wb=wallv0::blocks(ww,ww,hard);bg=background(p);current=solve(p,bg,pairs,&p,dt,&wb);record(current,ww,time>=tmax-1e-16);
  if(rank==0){history<<accepted<<','<<time<<','<<dt<<','<<maxres<<','<<(hard?-1e-12:min_wall)<<','<<sg<<','<<second.wall_constraints.size()<<','<<second.constraints.size()-second.wall_constraints.size()<<','<<rejected<<','<<ownership_changes<<'\n';if(accepted%100==0){history.flush();stages.flush();std::cout<<"accepted="<<accepted<<" time="<<time<<std::endl;}}
 }
 if(accepted>=cap&&time<tmax-1e-16)reason="STEP_CAP_NOT_COMPLETE";record(current,ww,true);
 if(rank==0){std::ofstream s("RUN_STATE.json");s<<std::setprecision(17)<<"{\"status\":\"TERMINAL\",\"reason\":\""<<reason<<"\",\"accepted_steps\":"<<accepted<<",\"time_s\":"<<time<<",\"target_time_s\":"<<tmax<<",\"max_residual\":"<<maxres<<",\"min_wall_gap_m\":"<<min_wall<<",\"min_pair_gap_m\":"<<min_pair<<",\"wall_constraint_activations\":"<<wall_activations<<",\"wall_active_observations\":"<<wall_active_observations<<",\"accepted_wall_overlap_count\":"<<wall_overlap<<",\"accepted_pair_overlap_count\":"<<pair_overlap<<",\"rejected_trials\":"<<rejected<<",\"invalid_query_count\":"<<invalid<<",\"nonfinite_count\":"<<nonfinite<<",\"lost_atoms\":0,\"owner_changes\":"<<ownership_changes<<",\"wall_limited_decisions\":"<<wall_limited<<",\"mpi_ranks\":"<<nr<<",\"kokkos_compatibility\":"<<(gpu?"true":"false")<<",\"hydrodynamics_enabled\":"<<(hydro?"true":"false")<<",\"hard_wall_enabled\":"<<(hard?"true":"false")<<",\"known_penetrating_control\":"<<(control?"true":"false")<<",\"RMBW_runtime_dependency\":false,\"physical_integrator\":\"RK2; unified wall/pair active-set; continuous swept acceptance\"}\n";std::cout<<"TERMINAL "<<reason<<" time="<<time<<std::endl;}
 delete l;l=nullptr;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();
 }catch(std::exception const&e){std::cerr<<"WALL_ENGINE_FAILED "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}
 MPI_Finalize();return 0;
}

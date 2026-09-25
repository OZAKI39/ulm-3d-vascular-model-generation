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
#include "adaptive_injection.hpp"
#include "workflow_geometry.hpp"
#include "lifecycle.hpp"
#include "kinematic_wall_constraint.hpp"
#include <functional>
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
static std::string number(double v){std::ostringstream o;o<<std::setprecision(17)<<v;return o.str();}
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,nr;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nr);LAMMPS*l=nullptr;
 try{
  if(argc!=2)throw std::runtime_error("USAGE workflow_lmp CONFIG");std::ifstream cfgfile(argv[1]);std::map<std::string,std::string>cfg;std::string k,v;while(cfgfile>>k>>v)cfg[k]=v;
  auto text=[&](std::string key){return cfg.at(key);};auto val=[&](std::string key){return std::stod(text(key));};
  std::string root=text("stage");double tmax=val("max_time"),dtmax=val("dt_max"),dtmin=val("dt_min"),margin=val("wall_margin"),size_lower=val("size_lower"),size_upper=val("size_upper");
  int stepscap=int(val("steps_cap")),maxretry=int(val("max_retry"));
  frozen::FlowField field(root+"/fields/FROZEN_FLOW_FIELD_V0.h5");FlowGradientSampler sampler(field);workflow::Geometry geometry(root);
  workflow::SurfacePatchTopology topology(root+"/geometry/WALL_ONLY.stl",root+"/geometry/WALL_TRIANGLE_REGIONS.txt");
  if(topology.nonmanifold_edges||topology.orientation_conflicts)throw std::runtime_error("STOP_WALL_SURFACE_TOPOLOGY_INVALID");
  double tie_factor=cfg.count("tie_factor")?val("tie_factor"):1,smooth_deg=cfg.count("smooth_dihedral_threshold_deg")?val("smooth_dihedral_threshold_deg"):15;
  int contact_step=0,contact_stage=0,contact_attempt=0;double contact_time=0,contact_dt=0;

  workflow::PositionSampler position(root,field,geometry,val("injection_clearance"),val("pair_clearance"),int(val("position_attempts")),text("position_mode")!="VALID_AREA_UNIFORM");
  workflow::NumberFluxProvider flux(val("authoritative_Q"),val("number_flux"),text("flux_mode")=="TEST_ONLY_DIRECT_NUMBER_FLUX");
  workflow::Controller controller;controller.source_basis=text("control_basis")=="SOURCE_POPULATION";controller.id_offset=0;
  std::vector<std::array<double,2>>size_stream;std::ifstream sizes(root+"/inputs/size_stream.txt");std::string line;while(std::getline(sizes,line)){if(line.empty()||line[0]=='#')continue;std::istringstream row(line);std::array<double,2>s;row>>s[0]>>s[1];size_stream.push_back(s);}
  std::vector<std::string>args={"workflow_lmp","-screen","none","-log","LAMMPS.log"};std::vector<char*>av;for(auto&s:args)av.push_back(s.data());l=new LAMMPS(av.size(),av.data(),MPI_COMM_WORLD);(*l->modify->fix_map)["sonovue/rigid/storage"]=&factory;
  auto command=[&](std::string c){l->input->one(c);};
  command("units si");command("atom_style sphere");command("boundary f f f");command("newton off");command("atom_modify map array");
  std::string box="region box block";for(int i=0;i<3;i++)box+=" "+number(field.origin[i]-20e-6)+" "+number(field.origin[i]+field.dx*field.dims[i]+20e-6);box+=" units box";command(box);command("create_box 1 box");
  command("pair_style zero 8e-6 full");command("pair_coeff * *");command("neighbor 4e-7 bin");command("neigh_modify every 1 delay 0 check no");command("comm_modify vel yes");command("fix rigid all sonovue/rigid/storage");command("thermo 1000000");command("run 0 post no");
  auto remove=[&](std::vector<long long>const&ids){workflow::remove_atoms(l,ids);};
  auto insert=[&](std::vector<workflow::Admission>const&admissions){
   int n=rank==0?admissions.size():0;MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);std::vector<double>data(n*5);if(!rank)for(int i=0;i<n;i++){data[5*i]=admissions[i].candidate.id;data[5*i+1]=admissions[i].candidate.radius;for(int j=0;j<3;j++)data[5*i+2+j]=admissions[i].x[j];}MPI_Bcast(data.data(),data.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
   std::vector<long long>out;std::vector<int>ids(n),types(n,1);std::vector<double>x(n*3);for(int i=0;i<n;i++){ids[i]=data[5*i];out.push_back(ids[i]);for(int j=0;j<3;j++)x[3*i+j]=data[5*i+2+j];}
   if(n){int made=lammps_create_atoms(l,n,ids.data(),types.data(),x.data(),nullptr,nullptr,0);if(made!=n)throw std::runtime_error("LAMMPS_INSERTION_COUNT_MISMATCH");std::map<int,double>radius;for(int i=0;i<n;i++)radius[ids[i]]=data[5*i+1];auto*a=l->atom;for(int i=0;i<a->nlocal;i++){auto it=radius.find(a->tag[i]);if(it!=radius.end()){a->radius[i]=it->second;a->rmass[i]=1000*4*std::acos(-1.)*std::pow(it->second,3)/3;for(int j=0;j<3;j++)a->omega[i][j]=0;}}command("run 0 post no");}
   return out;
  };
  std::ofstream traj,events,history,stages,flow,exitlog,nearlog,retrylog,pairlog,walllog,projectionlog,surfacelog;
  if(!rank){
   surfacelog.open("WALL_SURFACE_CONTACT_EVENTS.csv");
   surfacelog<<"time,step,stage,attempt,accepted,particle_id,radius,x,y,z,gap,raw_vx,raw_vy,raw_vz,used_vx,used_vy,used_vz,candidate_triangle_count,surface_cluster_count,constraint_mode,surface_mode,primary_surface_cluster,active_constraint_count,minimum_dihedral,maximum_dihedral,pseudonormal_x,pseudonormal_y,pseudonormal_z,dt,query_x,query_y,query_z,first_contact_fraction,candidate_ids,cluster_ids,active_indices,normal_set,safe,raw_safe,disclaimer\n";
   walllog.open("wall_constraint_timeseries.csv");projectionlog.open("GEOMETRIC_PROJECTION_EVENTS.csv");
   walllog<<"time,step,stage,particle_id,gap,gap_over_radius,nx,ny,nz,vx_raw,vy_raw,vz_raw,vx_used,vy_used,vz_used,vn_raw,vn_used,tangential_speed_raw,tangential_speed_used,removed_normal_speed,constraint_status,dt,wall_constraint_active,raw_segment_safe,omega_raw_x,omega_raw_y,omega_raw_z,omega_used_x,omega_used_y,omega_used_z,x_eval,y_eval,z_eval,disclaimer\n";
   projectionlog<<"particle_id,time,step,stage,attempt,accepted,gap_before,gap_after,correction_distance,nx,ny,nz,before_x,before_y,before_z,after_x,after_y,after_z,reason,disclaimer\n";
   traj.open("TRAJECTORIES.csv");events.open("injection_events.csv");history.open("SOLVER_HISTORY.csv");stages.open("INTEGRATION_STAGES.csv");flow.open("flux_timeseries.csv");exitlog.open("PORT_EVENTS.csv");nearlog.open("NEAR_WALL_EVENTS.csv");retrylog.open("RETRY_HISTORY.csv");pairlog.open("PAIR_HYDRODYNAMIC_EVENTS.csv");
   for(auto*f:{&traj,&events,&history,&stages,&flow,&exitlog,&nearlog,&retrylog,&pairlog,&walllog,&projectionlog,&surfacelog})*f<<std::setprecision(17);
   traj<<"step,time_s,particle_id,x_m,y_m,z_m,radius_m,vx_m_s,vy_m_s,vz_m_s,omega_x,omega_y,omega_z,wall_distance_m,surface_gap_m,wall_region,owner_rank,state,disclaimer\n";
   events<<"step,time_s,particle_id,event,radius_m,source_u,born_s,x_m,y_m,z_m,disclaimer\n";
   history<<"step,time_s,dt_s,active_before_removal,neighbor_pairs,hydrodynamic_pairs,residual_start,residual_mid,min_swept_gap_m,retries,crossrank_pairs,lammps_natoms,disclaimer\n";
   stages<<"step,stage,time_s,dt_s,particle_id,x_m,y_m,z_m,radius_m,vx_m_s,vy_m_s,vz_m_s,omega_x,omega_y,omega_z,disclaimer\n";
   flow<<"step,time,dt_s,authoritative_Q_in,target_number_flux,target_cumulative,admitted_cumulative,pending_count,active_count,outlet_0_cumulative,outlet_1_cumulative,outlet_2_cumulative,inlet_backflow_cumulative,size_rejected_cumulative,flux_deficit,N_source_drawn,N_terminal_failure,control_basis,disclaimer\n";
   exitlog<<"step,time_start_s,time_end_s,particle_id,event,port,fraction,start_x,start_y,start_z,end_x,end_y,end_z,cross_x,cross_y,cross_z,radius_m,disclaimer\n";
   nearlog<<"step,time_s,particle_id,event,gap_over_radius,wall_region,residence_s,disclaimer\n";
   retrylog<<"step,time_s,dt_attempt_s,trial_target,inserted_count,rolled_back_count,reason,disclaimer\n";
   pairlog<<"step,stage,time_s,id_i,id_j,gap_m,normal_resistance,shear_resistance,pump_resistance,force_x,force_y,force_z,constraint,twist,disclaimer\n";
  }
  auto write_stage=[&](int step,int stage,double time,double dt,const std::vector<Particle>&p,const Solution&q){if(rank)return;for(size_t i=0;i<p.size();i++){stages<<step<<','<<stage<<','<<time<<','<<dt<<','<<p[i].id;for(double x:p[i].x)stages<<','<<x;stages<<','<<p[i].a;for(int j=0;j<6;j++)stages<<','<<q.q[6*i+j];stages<<','<<workflow::disclaimer<<'\n';}for(auto&e:q.events){pairlog<<step<<','<<stage<<','<<time<<','<<p[e.i].id<<','<<p[e.j].id<<','<<e.gap<<','<<e.c.sq<<','<<e.c.sh<<','<<e.c.pu;for(double f:e.force)pairlog<<','<<f;pairlog<<','<<e.constraint<<",PENDING_NOT_IMPLEMENTED,"<<workflow::disclaimer<<'\n';}};
  auto write_wall=[&](int step,int stage,double t,double dt,const std::vector<Particle>&p,const Solution&raw,const Solution&used,const std::vector<workflow::WallConstraintResult>&wr){if(rank)return;for(size_t i=0;i<p.size();i++){auto const&r=wr[i];walllog<<t<<','<<step<<','<<stage<<','<<p[i].id<<','<<r.gap<<','<<r.gap/p[i].a;for(auto v:{r.normal,r.raw_velocity,r.constrained_velocity})for(double x:v)walllog<<','<<x;double vn=dot(r.constrained_velocity,r.normal);walllog<<','<<r.raw_normal_velocity<<','<<vn<<','<<norm(sub(r.raw_velocity,mul(r.normal,r.raw_normal_velocity)))<<','<<norm(sub(r.constrained_velocity,mul(r.normal,vn)))<<','<<r.removed_normal_velocity<<','<<r.status<<','<<dt<<','<<r.active<<','<<r.raw_segment_safe;for(int j=3;j<6;j++)walllog<<','<<raw.q[6*i+j];for(int j=3;j<6;j++)walllog<<','<<used.q[6*i+j];for(double x:p[i].x)walllog<<','<<x;walllog<<','<<workflow::disclaimer<<'\n';}};
  auto write_projection=[&](int step,int stage,int attempt,bool accepted,double t,const std::vector<Particle>&p,const std::vector<workflow::WallConstraintResult>&wr){if(rank)return;for(size_t i=0;i<wr.size();i++){auto const&r=wr[i];if(!r.position_projection_used)continue;projectionlog<<p[i].id<<','<<t<<','<<step<<','<<stage<<','<<attempt<<','<<accepted<<','<<r.gap_before<<','<<r.gap_after<<','<<r.position_projection_distance;for(auto v:{r.correction_normal,r.uncorrected_endpoint,r.endpoint})for(double x:v)projectionlog<<','<<x;projectionlog<<",BOUNDED_GEOMETRIC_OFFSET_PROJECTION,"<<workflow::disclaimer<<'\n';}};
  auto contact_failure=[&](std::string why,const Particle&p,Vec base,const workflow::WallConstraintResult&r,double h){if(rank)return;std::ofstream out("WALL_CONTACT_FAILURE.json");out<<std::setprecision(17)<<"{\"reason\":\""<<why<<"\",\"time\":"<<contact_time<<",\"step\":"<<contact_step<<",\"RK2_stage\":"<<contact_stage<<",\"attempt\":"<<contact_attempt<<",\"dt\":"<<contact_dt<<",\"stage_dt\":"<<h<<",\"particle_id\":"<<p.id<<",\"radius\":"<<p.a<<",\"position\":";workflow::write_vec(out,p.x);out<<",\"base\":";workflow::write_vec(out,base);out<<",\"raw_velocity\":";workflow::write_vec(out,r.raw_velocity);out<<",\"used_velocity\":";workflow::write_vec(out,r.constrained_velocity);out<<",\"endpoint\":";workflow::write_vec(out,r.endpoint);out<<",\"predicted_raw_swept_gap\":"<<geometry.wall.segment_distance_exact(base,add(base,mul(r.raw_velocity,h)))-p.a<<",\"used_swept_gap\":"<<geometry.wall.segment_distance_exact(base,r.endpoint)-p.a<<",\"tie_factor\":"<<tie_factor<<",\"smooth_dihedral_threshold_deg\":"<<smooth_deg<<",\"contact\":";topology.write_json(out,r.surface);out<<",\"endpoint_contact\":";if(r.endpoint_surface.candidates.empty())out<<"null";else topology.write_json(out,r.endpoint_surface);out<<",\"active_normal_indices\":[";for(size_t k=0;k<r.projection.active.size();k++){if(k)out<<',';out<<r.projection.active[k];}out<<"]}\n";out.flush();};
  auto write_surface=[&](int step,int stage,int attempt,bool accepted,double t,double h,const std::vector<Particle>&p,const std::vector<workflow::WallConstraintResult>&wr){if(rank)return;for(size_t i=0;i<wr.size();i++){auto const&r=wr[i];auto const&q=r.surface;surfacelog<<t<<','<<step<<','<<stage<<','<<attempt<<','<<accepted<<','<<p[i].id<<','<<p[i].a;for(double x:p[i].x)surfacelog<<','<<x;surfacelog<<','<<r.gap;for(auto v:{r.raw_velocity,r.constrained_velocity})for(double x:v)surfacelog<<','<<x;surfacelog<<','<<q.candidates.size()<<','<<q.clusters.size()<<','<<r.mode<<','<<q.status<<','<<(q.clusters.empty()?-1:q.clusters[0].id)<<','<<r.projection.active.size()<<','<<q.min_dihedral<<','<<q.max_dihedral;for(double x:r.normal)surfacelog<<','<<x;surfacelog<<','<<h;for(double x:r.contact_point)surfacelog<<','<<x;surfacelog<<','<<r.first_contact_fraction<<',';for(size_t j=0;j<q.candidates.size();j++){if(j)surfacelog<<'|';surfacelog<<q.candidates[j].id;}surfacelog<<',';for(size_t j=0;j<q.clusters.size();j++){if(j)surfacelog<<'|';surfacelog<<q.clusters[j].id;}surfacelog<<',';for(size_t j=0;j<r.projection.active.size();j++){if(j)surfacelog<<'|';surfacelog<<r.projection.active[j];}surfacelog<<',';for(size_t j=0;j<q.clusters.size();j++){if(j)surfacelog<<'|';for(int k=0;k<3;k++){if(k)surfacelog<<';';surfacelog<<q.clusters[j].normal[k];}}surfacelog<<','<<r.safe<<','<<r.raw_segment_safe<<','<<workflow::disclaimer<<'\n';}};
  auto constrain=[&](const std::vector<Particle>&eval,const std::vector<Particle>&base,const Solution&raw,Solution&used,std::vector<Particle>&end,std::vector<workflow::WallConstraintResult>&wr,double h,bool offset){used=raw;end=base;wr.clear();bool safe=true;for(size_t i=0;i<eval.size();i++){Vec v{raw.q[6*i],raw.q[6*i+1],raw.q[6*i+2]};workflow::WallConstraintResult r;try{r=workflow::constrain_surface_wall({eval[i].x,base[i].x,v,eval[i].a,margin,h,offset},geometry.wall,topology,tie_factor,smooth_deg);}catch(std::exception const&e){r.raw_velocity=v;r.constrained_velocity=v;r.endpoint=add(base[i].x,mul(v,h));r.surface=topology.query(eval[i].x,tie_factor,smooth_deg);contact_failure(e.what(),eval[i],base[i].x,r,h);throw;}if(r.status=="FAIL_NONMANIFOLD"||r.status.find("STOP_")==0){contact_failure(r.status,eval[i],base[i].x,r,h);throw std::runtime_error(r.status);}wr.push_back(r);for(int j=0;j<3;j++)used.q[6*i+j]=r.constrained_velocity[j];for(int j=3;j<6;j++)if(used.q[6*i+j]!=raw.q[6*i+j])throw std::runtime_error("FAIL_ROTATION_MODIFIED");end[i].x=r.endpoint;safe=safe&&r.safe;}return safe;};
  auto background=[&](const std::vector<Particle>&p){std::vector<Background>bg;for(auto&x:p)bg.push_back(sampler.query(x.x));return bg;};
  long long exits[4]={0,0,0,0},terminal=0,crossrank=0,total_inserted=0,total_removed=0,hydro_events=0;int accepted=0,retry_count=0;double time=0,min_swept=1e100,maxres=0,min_wall_gap=1e100;std::string reason="MAX_PHYSICAL_TIME";std::map<long long,double>near_start;
  auto series=[&](double dt,size_t active){if(rank)return;flow<<accepted<<','<<time<<','<<dt<<','<<flux.flow.volume_flux(time)<<','<<flux.number_flux(time)<<','<<controller.target<<','<<controller.admitted<<','<<controller.pending.size()<<','<<active<<','<<exits[1]<<','<<exits[2]<<','<<exits[3]<<','<<exits[0]<<','<<controller.size_rejected<<','<<controller.flux_debt()<<','<<controller.source_drawn<<','<<terminal<<','<<text("control_basis")<<','<<workflow::disclaimer<<'\n';flow.flush();};
  series(0,0);
  auto record=[&](int step,double t,std::vector<Particle>const&p,Solution const&q,std::vector<workflow::Crossing>const&crossings){if(rank)return;for(size_t i=0;i<p.size();i++){
   double wd=geometry.wall.distance(p[i].x),gap=wd-p[i].a;int region=geometry.region(p[i].x);min_wall_gap=std::min(min_wall_gap,gap);bool exited=crossings.size()&&crossings[i].port>=0;
   traj<<step<<','<<t<<','<<p[i].id;for(double x:p[i].x)traj<<','<<x;traj<<','<<p[i].a;for(int j=0;j<6;j++)traj<<','<<q.q[6*i+j];traj<<','<<wd<<','<<gap<<','<<region<<','<<p[i].owner<<','<<(exited?"EXITED":"ACTIVE")<<','<<workflow::disclaimer<<'\n';
   bool near=region==0&&gap/p[i].a<=val("nearwall_threshold");auto it=near_start.find(p[i].id);
   if(near&&it==near_start.end()){near_start[p[i].id]=t;nearlog<<step<<','<<t<<','<<p[i].id<<",ENTER,"<<gap/p[i].a<<','<<region<<",0,"<<workflow::disclaimer<<'\n';}
   if((!near||exited)&&it!=near_start.end()){nearlog<<step<<','<<t<<','<<p[i].id<<",LEAVE,"<<gap/p[i].a<<','<<region<<','<<t-it->second<<','<<workflow::disclaimer<<'\n';near_start.erase(it);}
  }};
  while(time<tmax-1e-15&&accepted<stepscap){
   auto original=particles(l,rank);double dt=std::min(dtmax,tmax-time);
   // The previous accepted population only proposes a dt. Injection is repeated
   // transactionally if new interactions tighten the displacement/gap bound.
   std::vector<Particle> p,mid,next;std::vector<std::pair<int,int>>pairs;Solution first,second,first_used,second_used;std::vector<workflow::WallConstraintResult>wallfirst,wallsecond;std::vector<workflow::Crossing>crossings;
   workflow::Controller trial;std::vector<workflow::InjectionEvent>injection_events;std::vector<long long>newids;bool ok=false;int retries=0;bool wall_retry=false;double sg=1e100;
   for(;retries<=maxretry&&dt>=dtmin;retries++){
    wallfirst.clear();wallsecond.clear();
    std::vector<workflow::Admission> admissions;injection_events.clear();
    if(!rank){trial=controller;auto draw=[&](long long i){return size_stream.at(i);};auto feasible=[&](workflow::Candidate const&c){if(c.radius>size_upper)return false;if(c.radius>size_lower)throw std::runtime_error("SIZE_CERTIFICATION_UNCERTAIN_BAND");return true;};auto place=[&](auto const&c,Vec&x,auto&r,auto const&fresh){return position.place(c,x,r,fresh,original);};admissions=trial.advance(time,dt,flux.number_flux(time),draw,feasible,place,injection_events,int(val("draw_budget")));}
    newids=insert(admissions);p=particles(l,rank);pairs=neighbors(l,p,crossrank);
    std::string rejection;double nextdt=dt*.5;
    try{
     auto bg=background(p);auto current=solve(p,bg,pairs);double allowed=dt;
     for(size_t i=0;i<p.size();i++){double speed=norm({current.q[6*i],current.q[6*i+1],current.q[6*i+2]});if(speed>0)allowed=std::min(allowed,.25*field.dx/speed);}
     for(auto&e:current.events){Vec rel{};for(int j=0;j<3;j++)rel[j]=current.q[6*e.i+j]-current.q[6*e.j+j];double re=p[e.i].a*p[e.j].a/(p[e.i].a+p[e.j].a);allowed=std::min(allowed,.4*std::max(e.gap,.001*re)/std::max(std::abs(dot(rel,e.n)),1e-30));}
     if(allowed<dt*(1-1e-12)){rejection="ADAPTIVE_DISPLACEMENT_OR_GAP";nextdt=allowed;}
     else{
      first=solve(p,bg,pairs,&p,.5*dt);
      contact_step=accepted;contact_stage=0;contact_time=time;contact_attempt=retries;contact_dt=dt;
      bool safe=constrain(p,p,first,first_used,mid,wallfirst,.5*dt,wall_retry);
      if(!safe)rejection="CURVED_GEOMETRY_MIDPOINT_RETRY";
      else if(swept(p,mid,pairs)<-1e-12)rejection="PAIR_SWEPT_MIDPOINT_REJECT";
      else{
       second=solve(mid,background(mid),pairs,&p,dt);
       contact_stage=1;contact_time=time+.5*dt;
       safe=constrain(mid,p,second,second_used,next,wallsecond,dt,wall_retry);
       sg=swept(p,next,pairs);if(!safe)rejection="CURVED_GEOMETRY_FINAL_RETRY";else if(sg< -1e-12)rejection="PAIR_SWEPT_REJECT";
       else{
        crossings.clear();for(size_t i=0;i<p.size();i++)crossings.push_back(geometry.crossing(p[i].x,next[i].x));
        // Exit proposals require no query outside the cap. All remaining active
        // centers still require the unchanged eight-corner native sampler.
        for(size_t i=0;i<p.size();i++)if(crossings[i].port<0)sampler.query(next[i].x);
        ok=true;
       }
      }
     }
    }catch(std::runtime_error const&e){std::string err=e.what();if(err=="INVALID_FLOW_QUERY"||err=="INVALID_GRADIENT_CORNER")rejection="FLOW_QUERY_REJECT";else throw;}
    write_projection(accepted,0,retries,ok,time,p,wallfirst);write_projection(accepted,1,retries,ok,time+.5*dt,mid,wallsecond);
    write_surface(accepted,0,retries,ok,time,.5*dt,p,wallfirst);write_surface(accepted,1,retries,ok,time+.5*dt,dt,mid,wallsecond);
    if(ok)break;
    if(rejection.find("CURVED_GEOMETRY")!=std::string::npos)wall_retry=true;
    remove(newids);auto restored=particles(l,rank);if(restored.size()!=original.size())throw std::runtime_error("ROLLBACK_PARTICLE_COUNT");for(size_t i=0;i<restored.size();i++)if(restored[i].id!=original[i].id||restored[i].x!=original[i].x)throw std::runtime_error("ROLLBACK_STATE_DRIFT");
    if(!rank)retrylog<<accepted<<','<<time<<','<<dt<<','<<trial.target<<','<<newids.size()<<','<<newids.size()<<','<<rejection<<','<<workflow::disclaimer<<'\n';retry_count++;
    reason=rejection=="FLOW_QUERY_REJECT"?"FLOW_QUERY_FAILED":(rejection.find("PAIR")!=std::string::npos?"PAIR_SAFETY_FAILED":"CURVED_GEOMETRY_TANGENTIAL_STALL");dt=nextdt;
   }
   if(!ok){if(reason=="CURVED_GEOMETRY_TANGENTIAL_STALL"){bool written=false;for(int stage=0;stage<2&&!written;stage++){auto const&wr=stage?wallsecond:wallfirst;auto const&eval=stage?mid:p;for(size_t i=0;i<wr.size();i++)if(!wr[i].safe){contact_stage=stage;contact_failure(reason,eval[i],p[i].x,wr[i],stage?contact_dt:.5*contact_dt);written=true;break;}}}break;}
   reason="MAX_PHYSICAL_TIME";if(!rank){controller=trial;for(auto&e:injection_events){events<<accepted<<','<<e.time<<','<<e.id<<','<<e.kind<<','<<e.radius<<','<<e.u<<','<<e.born;for(double x:e.x)events<<','<<x;events<<','<<workflow::disclaimer<<'\n';}}
   total_inserted+=newids.size();write_stage(accepted,0,time,dt,p,first);write_stage(accepted,1,time+.5*dt,dt,mid,second);
   write_wall(accepted,0,time,dt,p,first,first_used,wallfirst);write_wall(accepted,1,time+.5*dt,dt,mid,second,second_used,wallsecond);
   // Admission starts at t_n and participates in both RK2 solves immediately.
   if(!newids.empty())record(accepted,time,p,first_used,{});
   proposal.clear();for(size_t i=0;i<p.size();i++){std::array<double,9>a;for(int j=0;j<3;j++){a[j]=next[i].x[j];a[j+3]=second_used.q[6*i+j];a[j+6]=second.q[6*i+j+3];}proposal[p[i].id]=a;}
   command("timestep "+number(dt));command("run 1 pre no post no");accepted++;double oldtime=time;time+=dt;
   auto actual=particles(l,rank);if(actual.size()!=next.size())throw std::runtime_error("FAIL_PARTICLE_ACCOUNTING_LAMMPS_LOSS");for(size_t i=0;i<actual.size();i++)if(actual[i].id!=next[i].id||norm(sub(actual[i].x,next[i].x))>1e-20)throw std::runtime_error("LAMMPS_STORAGE_RK2_MISMATCH");
   record(accepted,time,actual,second_used,crossings);std::vector<long long>departed;
   for(size_t i=0;i<actual.size();i++)if(crossings[i].port>=0){auto&c=crossings[i];exits[c.port]++;departed.push_back(actual[i].id);if(!rank){exitlog<<accepted<<','<<oldtime<<','<<time<<','<<p[i].id<<','<<(c.port?"EXITED_OUTLET":"EXITED_INLET_BACKFLOW")<<','<<geometry.ports[c.port].name<<','<<c.fraction;for(auto x:{p[i].x,next[i].x,c.x})for(double v:x)exitlog<<','<<v;exitlog<<','<<p[i].a<<','<<workflow::disclaimer<<'\n';}}
   remove(departed);total_removed+=departed.size();auto active=particles(l,rank);min_swept=std::min(min_swept,sg);maxres=std::max({maxres,first.residual,second.residual});hydro_events+=first.events.size()+second.events.size();
   if(!rank){if(controller.admitted!=static_cast<long long>(active.size())+exits[0]+exits[1]+exits[2]+exits[3]+terminal)throw std::runtime_error("FAIL_PARTICLE_ACCOUNTING");history<<accepted<<','<<time<<','<<dt<<','<<p.size()<<','<<pairs.size()<<','<<second.events.size()<<','<<first.residual<<','<<second.residual<<','<<sg<<','<<retries<<','<<crossrank<<','<<l->atom->natoms<<','<<workflow::disclaimer<<'\n';}
   series(dt,active.size());int capacity=!rank&&controller.capacity_exceeded(time,size_t(val("max_pending")),val("pending_max_age"));MPI_Bcast(&capacity,1,MPI_INT,0,MPI_COMM_WORLD);if(capacity){reason="INJECTION_CAPACITY_EXCEEDED";break;}
   if(!rank&&accepted%500==0){traj.flush();events.flush();stages.flush();history.flush();surfacelog.flush();std::cout<<"accepted="<<accepted<<" time="<<time<<" active="<<active.size()<<" admitted="<<controller.admitted<<" pending="<<controller.pending.size()<<std::endl;}
  }
  if(accepted>=stepscap&&time<tmax-1e-15)reason="STEP_CAP_NOT_COMPLETE";
  if(!rank){for(auto&n:near_start)nearlog<<accepted<<','<<time<<','<<n.first<<",RIGHT_CENSORED,0,0,"<<time-n.second<<','<<workflow::disclaimer<<'\n';
   std::ofstream out("RUN_STATE.json");out<<std::setprecision(17)<<"{\"status\":\"TERMINAL\",\"reason\":\""<<reason<<"\",\"time_s\":"<<time<<",\"accepted_steps\":"<<accepted<<",\"source_drawn\":"<<controller.source_drawn<<",\"size_rejected\":"<<controller.size_rejected<<",\"admitted\":"<<controller.admitted<<",\"pending\":"<<controller.pending.size()<<",\"active\":"<<l->atom->natoms<<",\"outlet_exits\":"<<exits[1]+exits[2]+exits[3]<<",\"inlet_backflow_exits\":"<<exits[0]<<",\"terminal_failures\":"<<terminal<<",\"max_residual\":"<<maxres<<",\"min_swept_gap_m\":"<<min_swept<<",\"min_wall_gap_m\":"<<min_wall_gap<<",\"hydrodynamic_pair_observations\":"<<hydro_events<<",\"retries\":"<<retry_count<<",\"committed_LAMMPS_insertions\":"<<total_inserted<<",\"committed_LAMMPS_removals\":"<<total_removed<<",\"mpi_ranks\":"<<nr<<",\"physics\":\"bulk + stable rigid bubble-bubble excess; twist pending\",\"wall_hydrodynamics\":false,\"neighbor_source\":\"LAMMPS pair_zero full list\",\"integrator\":\"stable rigid midpoint RK2 with LAMMPS storage and migration\",\"disclaimer\":\""<<workflow::disclaimer<<"\"}\n";
   std::ofstream pending("PENDING_FINAL.csv");pending<<"particle_id,radius_m,born_s,source_u,disclaimer\n"<<std::setprecision(17);for(auto&c:controller.pending)pending<<c.id<<','<<c.radius<<','<<c.born<<','<<c.u<<','<<workflow::disclaimer<<'\n';std::cout<<"TERMINAL "<<reason<<" time="<<time<<std::endl;
  }
  delete l;l=nullptr;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();MPI_Finalize();return 0;
 }catch(std::exception const&e){std::cerr<<"WORKFLOW_FAILED "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}

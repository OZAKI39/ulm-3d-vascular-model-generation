#define main unused_volume_harness_entry
#include "native_harness.cpp"
#undef main
#include "Simulation.h"
#include "all_fun.h"
#include "read_files.h"
#include "distribute.h"
#include "baf_ini.h"
#include "set_bc.h"
#include "dispatch_trace.h"

void allocate(SolutionStates& s,int n) {
 for(auto* state:{&s.old,&s.current,&s.intermediate}) {
  state->get_velocity().resize(4,n); state->get_acceleration().resize(4,n); state->get_displacement().resize(4,n);
 }
}
int main(int argc,char** argv) {
 MPI_Init(&argc,&argv);
 try {
  if(argc!=3)throw std::runtime_error("XML and output-prefix required");
  std::string prefix=argv[2];dispatch_trace::out.open(prefix+"_trace.csv");dispatch_trace::out<<"stage,face,function\n";
  Simulation sim;read_files_ns::read_files(&sim,argv[1]);distribute(&sim);
  auto& cm=sim.com_mod;cm.cEq=0;cm.cDmn=0;cm.tDof=4;cm.dof=4;cm.Bf.resize(3,cm.tnNo);
  auto& eq=cm.eq[0];eq.dof=4;eq.s=0;eq.e=3;eq.af=2./3;eq.am=5./6;eq.gam=2./3;
  SolutionStates sol;allocate(sol,cm.tnNo);Capture cap;eq.linear_algebra=&cap;
  std::ofstream props(prefix+"_properties.csv");props<<std::setprecision(17)<<"face,bType,weakDir,face_nodes,face_elements,general_dof,gx_size,gx_zero,flux_integral,pressure_changed,max_wall_u\n";
  for(auto& mesh:cm.msh) {
   fs::init_fs_msh(cm,mesh);
   for(auto& face:mesh.fa){fs::init_fs_face(cm,mesh,face);baf_ini_ns::face_ini(&sim,mesh,face,sol);}
  }
  for(auto& bc:eq.bc){auto& face=cm.msh[bc.iM].fa[bc.iFa];baf_ini_ns::bc_ini(cm,sim.cm_mod,bc,face,sol);}
  for(int a=0;a<cm.tnNo;a++)sol.current.get_velocity()(3,a)=7.;
  dispatch_trace::stage="actual_P1";
  set_bc::set_bc_dir(cm,sol);set_bc::set_bc_dir_w(cm,sol);set_bc::set_bc_neu(cm,sim.cm_mod,sol);
  for(auto& bc:eq.bc) {
   auto& face=cm.msh[bc.iM].fa[bc.iFa];int zeros=0;double flux=0.,pc=0.,wall=0.;Vector<double> gx(cm.tnNo);
   for(int a=0;a<bc.gx.size();a++){zeros+=bc.gx(a)==0.;gx(face.gN(a))=bc.gx(a);}
   if(bc.gx.size())flux=all_fun::integ(cm,sim.cm_mod,face,gx,sol,false);
   for(int a=0;a<face.nNo;a++){int id=face.gN(a);pc=std::max(pc,std::abs(sol.current.get_velocity()(3,id)-7.));for(int j=0;j<3;j++)wall=std::max(wall,std::abs(sol.current.get_velocity()(j,id)));}
   props<<face.name<<','<<bc.bType<<','<<bc.weakDir<<','<<face.nNo<<','<<face.nEl<<','<<bc.gm.dof<<','<<bc.gx.size()<<','<<zeros<<','<<flux<<','<<pc<<','<<wall<<'\n';
  }
  // Replay native parser-created objects on a curved TET10/TRI6 fixture.
  // No type/weak flag is assigned by this harness. Only geometry/data dimensions
  // are localized, after actual native parser + full target profile initialization.
  dispatch_trace::stage="mixed_P2_replay";
  std::ofstream mixed(prefix+"_mixed_checks.csv");mixed<<std::setprecision(17)<<"face,strong_all_six_error,pressure_error,unconstrained_error,neumann_R_scaled_error,neumann_K_max\n";
  for(const auto& sourcebc:eq.bc) {
   const auto& sourceface=cm.msh[sourcebc.iM].fa[sourcebc.iFa];
   Simulation mini;auto& c=mini.com_mod;c.nsd=3;c.nsymd=6;c.nEq=1;c.nMsh=1;c.tDof=4;c.dof=4;c.tnNo=10;c.cEq=0;c.cDmn=0;c.dt=.01;c.mvMsh=false;
   c.eq.resize(1);c.eq[0]=eq;c.eq[0].nBc=1;c.eq[0].bc={sourcebc};c.eq[0].linear_algebra=&cap;
   c.msh.resize(1);c.x.resize(3,10);c.Bf.resize(3,10);
   const double v[4][3]={{.001,0,0},{0,.001,0},{0,0,.001},{0,0,0}};
   for(int a=0;a<4;a++)for(int j=0;j<3;j++)c.x(j,a)=v[a][j];
   const int edges[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
   for(int k=0;k<6;k++)for(int j=0;j<3;j++)c.x(j,k+4)=(v[edges[k][0]][j]+v[edges[k][1]][j])*.5;
   c.x(2,4)=.00001; // positively oriented curved geometry, independent of dispatch
   auto& m=c.msh[0];m.eNoN=10;m.nFs=2;m.nEl=1;m.nFa=1;m.lShl=false;m.lFib=false;m.IEN.resize(10,1);
   for(int a=0;a<10;a++)m.IEN(a,0)=a;
   nn::select_ele(c,m);fs::init_fs_msh(c,m);m.fa.resize(1);auto& f=m.fa[0];f.name=sourceface.name;f.iM=0;f.eNoN=6;f.nNo=6;f.nEl=1;f.nG=7;f.eType=consts::ElementType::TRI6;
   const int ids[6]={0,1,3,4,8,7};f.IEN.resize(6,1);f.gN.resize(6);f.gE.resize(1);f.gE(0)=0;f.nV.resize(3,6);
   for(int a=0;a<6;a++){f.IEN(a,0)=ids[a];f.gN(a)=ids[a];f.nV(2,a)=-1.;}
   f.w.resize(7);f.xi.resize(2,7);f.N.resize(6,7);f.Nx.resize(2,6,7);nn::get_gip(nullptr,f);for(int g=0;g<7;g++)nn::get_gnn(nullptr,g,f);
   fs::init_fs_face(c,m,f);
   auto& bc=c.eq[0].bc[0];bc.iM=0;bc.iFa=0;
   if(sourcebc.gx.size()){bc.gx.resize(6);for(int a=0;a<6;a++)bc.gx(a)=sourcebc.gx(a%sourcebc.gx.size());}
   if(sourcebc.gm.dof>0 && sourcebc.gm.nTP>0){bc.gm.d.resize(sourcebc.gm.dof,6,sourcebc.gm.nTP);for(int t=0;t<bc.gm.nTP;t++)for(int a=0;a<6;a++)for(int j=0;j<bc.gm.dof;j++)bc.gm.d(j,a,t)=sourcebc.gm.d(j,a%sourceface.nNo,t);}
   SolutionStates ss;allocate(ss,10);for(int a=0;a<10;a++)for(int j=0;j<4;j++)ss.current.get_velocity()(j,a)=7.;
   set_bc::set_bc_dir(c,ss);set_bc::set_bc_dir_w(c,ss);set_bc::set_bc_neu(c,mini.cm_mod,ss);
   double pressure_error=0.;for(int a=0;a<10;a++)pressure_error=std::max(pressure_error,std::abs(ss.current.get_velocity()(3,a)-7.));
   if(pressure_error!=0.)throw std::runtime_error("Boundary modified pressure slots");
   double strong_error=0.,untouched=0.,neu_error=0.,kmax=0.;
   bool strong=utils::btest(bc.bType,consts::iBC_Dir)&&!bc.weakDir;
   for(int a=0;a<10;a++) {
    int onface=-1;for(int k=0;k<6;k++)if(ids[k]==a)onface=k;
    for(int j=0;j<3;j++) {
     double value=7.;
     if(strong&&onface>=0) value=utils::btest(bc.bType,consts::iBC_gen)?bc.gm.d(j,onface,0):bc.g*bc.gx(onface)*f.nV(j,onface);
     double err=std::abs(ss.current.get_velocity()(j,a)-value);
     if(strong&&onface>=0)strong_error=std::max(strong_error,err);else untouched=std::max(untouched,err);
    }
   }
   if(utils::btest(bc.bType,consts::iBC_Neu)) {
    Array<double> ref(4,6);
    for(int g=0;g<7;g++) {
     double tangent[3][2]={};
     for(int a=0;a<6;a++)for(int j=0;j<3;j++)for(int k=0;k<2;k++)tangent[j][k]+=f.Nx(k,a,g)*c.x(j,ids[a]);
     double normal[3]={tangent[1][0]*tangent[2][1]-tangent[2][0]*tangent[1][1],tangent[2][0]*tangent[0][1]-tangent[0][0]*tangent[2][1],tangent[0][0]*tangent[1][1]-tangent[1][0]*tangent[0][1]};
     double dot=0.;for(int j=0;j<3;j++)dot+=normal[j]*(c.x(j,0)-c.x(j,2));
     if(dot<0.)for(int j=0;j<3;j++)normal[j]*=-1.;
     double h=0.;for(int a=0;a<6;a++)h-=f.N(a,g)*bc.g*bc.gx(a);
     for(int a=0;a<6;a++)for(int j=0;j<3;j++)ref(j,a)-=f.w(g)*f.N(a,g)*h*normal[j];
    }
    double scale=1e-30;for(int a=0;a<6;a++)for(int j=0;j<4;j++){scale=std::max(scale,std::abs(ref(j,a)));neu_error=std::max(neu_error,std::abs(ref(j,a)-cap.R(j,a)));}
    neu_error/=scale;
    for(int a=0;a<6;a++)for(int b=0;b<6;b++)for(int j=0;j<16;j++)kmax=std::max(kmax,std::abs(cap.K(j,a,b)));
   }
   mixed<<f.name<<','<<strong_error<<','<<pressure_error<<','<<untouched<<','<<neu_error<<','<<kmax<<'\n';
   if(strong_error>1e-11||untouched!=0.||neu_error>1e-11||kmax!=0.)throw std::runtime_error("Mixed boundary numerical regression failed");
   c.eq[0].linear_algebra=nullptr;
  }
  eq.linear_algebra=nullptr;
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Finalize();return 1;}
 MPI_Finalize();return 0;
}

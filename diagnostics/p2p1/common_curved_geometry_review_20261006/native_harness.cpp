// Actual native construct_fluid and retained element kernels. Never solve/step.
#include "nn.h"
#include "fs.h"
#include "fluid.h"
#include "LinearAlgebra.h"
#include "trace.h"
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <stdexcept>

struct Capture: LinearAlgebra {
  Array<double> R; Array3<double> K;
  void alloc(ComMod&,eqType&) override {}
  void assemble(ComMod&,int,const Vector<int>&,const Array3<double>& k,const Array<double>& r) override {R=r;K=k;}
  void check_options(consts::PreconditionerType,consts::LinearAlgebraType) override {}
  void initialize(ComMod&,eqType&) override {}
  void set_assembly(consts::LinearAlgebraType) override {}
  void set_preconditioner(consts::PreconditionerType) override {}
  void solve(ComMod&,eqType&,const Vector<int>&,const Vector<double>&) override {throw std::runtime_error("Linear solves forbidden");}
  void finalize() override {}
};
static void read_array(std::istream& in,Array<double>& a) {
  for(int n=0;n<a.ncols();++n)for(int j=0;j<a.nrows();++j)
    if(!(in>>a(j,n)))throw std::runtime_error("Incomplete input array");
}
static void write_capture(const std::string& path,const Capture& c,int n) {
  std::ofstream out(path);out<<std::setprecision(17);
  for(int i=0;i<4;i++)for(int a=0;a<n;a++)out<<c.R(i,a)<<' ';
  out<<'\n';
  for(int i=0;i<16;i++)for(int a=0;a<n;a++)for(int b=0;b<n;b++)out<<c.K(i,a,b)<<' ';
  out<<'\n';
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  try {
    if(argc<4)throw std::runtime_error("input output-prefix native|reference [reference-file]");
    std::ifstream in(argv[1]);int n;in>>n;if(n!=4 && n!=10)throw std::runtime_error("TET4 or TET10 only");
    ComMod cm;cm.nsd=3;cm.nsymd=6;cm.tDof=4;cm.dof=4;cm.cEq=0;cm.cDmn=0;cm.dt=.01;cm.nMsh=1;cm.tnNo=n;
    cm.urisFlag=false;cm.risFlag=false;cm.eq.resize(1);cm.msh.resize(1);cm.x.resize(3,n);cm.Bf.resize(3,n);
    read_array(in,cm.x);
    auto& eq=cm.eq[0];eq.phys=consts::EquationType::phys_fluid;eq.dof=4;eq.s=0;eq.e=3;eq.nDmn=1;
    eq.af=2./3;eq.am=5./6;eq.gam=2./3;eq.dmn.resize(1);
    auto& dm=eq.dmn[0];dm.Id=-1;dm.phys=consts::EquationType::phys_fluid;
    dm.prop[consts::PhysicalPropertyType::fluid_density]=1056.;
    dm.prop[consts::PhysicalPropertyType::brinkman_inverse_permeability]=0.;
    dm.prop[consts::PhysicalPropertyType::backflow_stab]=0.;
    dm.fluid_visc.viscType=consts::FluidViscosityModelType::viscType_Const;dm.fluid_visc.mu_i=.00345312;
    auto& mesh=cm.msh[0];mesh.nFs=n==10?2:1;mesh.nEl=1;mesh.eNoN=n;mesh.lShl=false;mesh.lFib=false;
    mesh.IEN.resize(n,1);for(int a=0;a<n;a++)mesh.IEN(a,0)=a;
    nn::select_ele(cm,mesh);fs::init_fs_msh(cm,mesh);
    SolutionStates state;state.intermediate.get_velocity().resize(4,n);state.intermediate.get_acceleration().resize(4,n);
    read_array(in,state.intermediate.get_velocity());read_array(in,state.intermediate.get_acceleration());read_array(in,cm.Bf);
    Capture cap;eq.linear_algebra=&cap;std::string prefix(argv[2]),mode(argv[3]);
    if(mode=="native") {
      validation_trace::points.open(prefix+"_q.txt");validation_trace::points<<std::setprecision(17);
      validation_trace::parameters.open(prefix+"_tau.txt");validation_trace::parameters<<std::setprecision(17);
      fluid::construct_fluid(cm,mesh,state);
    } else if(mode=="reference" && argc==5) {
      // Independent Python geometry and transformations; same unmodified native
      // fluid physics kernels and quadrature. No alternative solver.
      std::ifstream ref(argv[4]);int count;ref>>count;
      cap.R.resize(4,n);cap.K.resize(16,n,n);
      for(int q=0;q<count;q++) {
        int opt;double w;ref>>opt>>w;
        Vector<double> Nw(n),Np(4),uris(3);Array<double> K(3,3),Dw(3,n),Dp(3,4),Hw(6,n);
        for(int a=0;a<n;a++)ref>>Nw(a);for(int a=0;a<4;a++)ref>>Np(a);
        read_array(ref,Dw);read_array(ref,Dp);read_array(ref,Hw);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)ref>>K(i,j);
        if(!ref)throw std::runtime_error("Truncated common-reference input");
        if(opt==1)fluid::fluid_3d_m(cm,n==4,n,4,w,K,Nw,Np,Dw,Dp,Hw,
            state.intermediate.get_acceleration(),state.intermediate.get_velocity(),cm.Bf,cap.R,cap.K,0.,0.,uris);
        else fluid::fluid_3d_c(cm,n==4,n,4,w,K,Nw,Np,Dw,Dp,Hw,
            state.intermediate.get_acceleration(),state.intermediate.get_velocity(),cm.Bf,cap.R,cap.K,0.,0.,uris);
      }
    } else throw std::runtime_error("Unknown validation mode");
    write_capture(prefix+"_RK.txt",cap,n);eq.linear_algebra=nullptr;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Finalize();return 1;}
  MPI_Finalize();return 0;
}

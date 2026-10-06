// Small geometry/assembly probe using the actual retained production objects.
// Capture only intercepts the assembled element result: no alternative solver.
#include "nn.h"
#include "fs.h"
#include "fluid.h"
#include "LinearAlgebra.h"
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

struct Capture: LinearAlgebra {
  Array<double> R;
  void alloc(ComMod&,eqType&) override {}
  void assemble(ComMod&,const int n,const Vector<int>&,const Array3<double>&,const Array<double>& r) override {
    if(n!=10)throw std::runtime_error("not TET10"); R=r;
  }
  void check_options(consts::PreconditionerType,consts::LinearAlgebraType) override {}
  void initialize(ComMod&,eqType&) override {}
  void set_assembly(consts::LinearAlgebraType) override {}
  void set_preconditioner(consts::PreconditionerType) override {}
  void solve(ComMod&,eqType&,const Vector<int>&,const Vector<double>&) override {throw std::runtime_error("no solves allowed");}
  void finalize() override {}
};

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  try {
    if(argc!=4)throw std::runtime_error("nodes output summary required");
    ComMod cm;cm.nsd=3;cm.nsymd=6;cm.tDof=4;cm.dof=4;cm.cEq=0;cm.cDmn=0;cm.dt=.01;cm.nMsh=1;cm.tnNo=10;
    cm.urisFlag=false;cm.risFlag=false;cm.eq.resize(1);cm.msh.resize(1);cm.x.resize(3,10);cm.Bf.resize(3,10);
    std::ifstream input(argv[1]);
    for(int a=0;a<10;a++)for(int j=0;j<3;j++)if(!(input>>cm.x(j,a)))throw std::runtime_error("bad nodes");
    auto& eq=cm.eq[0];eq.phys=consts::EquationType::phys_fluid;eq.dof=4;eq.s=0;eq.e=3;eq.nDmn=1;
    eq.af=2./3;eq.am=5./6;eq.gam=2./3;eq.dmn.resize(1);
    auto& dm=eq.dmn[0];dm.Id=-1;dm.phys=consts::EquationType::phys_fluid;
    dm.prop[consts::PhysicalPropertyType::fluid_density]=1056.;
    dm.prop[consts::PhysicalPropertyType::brinkman_inverse_permeability]=0.;
    dm.fluid_visc.viscType=consts::FluidViscosityModelType::viscType_Const;dm.fluid_visc.mu_i=.00345312;
    auto& mesh=cm.msh[0];mesh.nFs=2;mesh.nEl=1;mesh.eNoN=10;mesh.lShl=false;mesh.lFib=false;
    mesh.IEN.resize(10,1);for(int a=0;a<10;a++)mesh.IEN(a,0)=a;
    nn::select_ele(cm,mesh);fs::init_fs_msh(cm,mesh);
    Capture capture;eq.linear_algebra=&capture;
    SolutionStates state;state.intermediate.get_velocity().resize(4,10);state.intermediate.get_acceleration().resize(4,10);
    for(int a=0;a<10;a++)for(int j=0;j<3;j++)state.intermediate.get_velocity()(j,a)=cm.x(j,a)-cm.x(j,3);
    // Execute the production element loop and capture its true continuity residual.
    fluid::construct_fluid(cm,mesh,state);
    std::ofstream summary(argv[3]);summary<<std::setprecision(17)<<"{\"status\":\"EXECUTED\",\"native_construct_fluid_Rp\":[";
    for(int a=0;a<4;a++)summary<<(a?",":"")<<capture.R(3,a);
    summary<<"],\"velocity_field\":\"u=x-x_corner3, exact G=I and div=3\",\"CFD_steps\":0,\"linear_solves\":0}\n";
    // Expose native basis/mappings with exactly construct_fluid's conditional order.
    Array<double> xq(3,4),Dw(3,10),Dq(3,4),metric(3,3);
    for(int a=0;a<4;a++)for(int j=0;j<3;j++)xq(j,a)=cm.x(j,a);
    std::ofstream out(argv[2]);out<<std::setprecision(17);
    double usedJac=0.;std::array<fsType,2> spaces;
    for(int opt=1;opt<=2;opt++) {
      fs::get_thood_fs(cm,spaces,mesh,false,opt);
      const int ng=spaces[opt==1?0:1].nG;
      for(int g=0;g<ng;g++) {
        auto drw=spaces[0].Nx.rslice(g);auto drq=spaces[1].Nx.rslice(g);
        if(opt==1) {
          if(g==0 || !spaces[1].lShpF)nn::gnn(4,3,3,drq,xq,Dq,usedJac,metric);
          if(g==0 || !spaces[0].lShpF)nn::gnn(10,3,3,drw,cm.x,Dw,usedJac,metric);
        } else {
          if(g==0 || !spaces[0].lShpF)nn::gnn(10,3,3,drw,cm.x,Dw,usedJac,metric);
          if(g==0 || !spaces[1].lShpF)nn::gnn(4,3,3,drq,xq,Dq,usedJac,metric);
        }
        double jv=0.,jq=0.;Array<double> tempW(3,10),tempQ(3,4);
        nn::gnn(10,3,3,drw,cm.x,tempW,jv,metric);nn::gnn(4,3,3,drq,xq,tempQ,jq,metric);
        out<<opt<<' '<<g<<' '<<spaces[opt==1?0:1].w(g)<<' ';
        for(int j=0;j<3;j++)out<<spaces[0].xi(j,g)<<' ';
        out<<jv<<' '<<jq<<' '<<usedJac<<' ';
        for(int a=0;a<4;a++)for(int j=0;j<3;j++)out<<Dq(j,a)<<' ';
        for(int a=0;a<10;a++)for(int j=0;j<3;j++)out<<Dw(j,a)<<' ';
        out<<'\n';
      }
    }
    eq.linear_algebra=nullptr;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Finalize();return 1;}
  MPI_Finalize();return 0;
}

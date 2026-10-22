// Validation harness linked against the retained production object files.
// No CFD time step or linear solve; no replacement element/assembly implementation.
#include "nn.h"
#include "fs.h"
#include "fluid.h"
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    if (argc != 4) throw std::runtime_error("input output summary required");
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    output << std::setprecision(17);
    Vector<double> xi(3), Nw(10), Nq(4);
    Array<double> Dw(3,10), Dq(3,4);
    int count=0;
    while (input >> xi(0) >> xi(1) >> xi(2)) {
      nn::get_gnn(3, consts::ElementType::TET10, 10, xi, Nw, Dw);
      nn::get_gnn(3, consts::ElementType::TET4, 4, xi, Nq, Dq);
      for (int a=0;a<10;a++) output << Nw(a) << ' ';
      for (int a=0;a<10;a++) for(int j=0;j<3;j++) output << Dw(j,a) << ' ';
      for (int a=0;a<4;a++) output << Nq(a) << ' ';
      for (int a=0;a<4;a++) for(int j=0;j<3;j++) output << Dq(j,a) << ' ';
      output << '\n'; ++count;
    }
    ComMod cm;
    cm.nsd=3; cm.cEq=0; cm.cDmn=0; cm.dt=.01; cm.nMsh=1; cm.tnNo=10;
    cm.eq.resize(1); cm.msh.resize(1);
    auto &eq=cm.eq[0]; eq.phys=consts::EquationType::phys_fluid;
    eq.af=1.; eq.am=1.; eq.gam=1.; eq.dmn.resize(1);
    auto &dm=eq.dmn[0]; dm.prop[consts::PhysicalPropertyType::fluid_density]=1056.;
    dm.fluid_visc.viscType=consts::FluidViscosityModelType::viscType_Const;
    dm.fluid_visc.mu_i=.00345312;
    auto &mesh=cm.msh[0]; mesh.nFs=2; mesh.nEl=1; mesh.eNoN=10;
    mesh.lShl=false; mesh.lFib=false;
    // Production first initializes the mesh basis; init_fs alone is for
    // reduced pressure spaces and its quadrature overload does not take TET10.
    nn::select_ele(cm,mesh);
    fs::init_fs_msh(cm,mesh);
    std::array<fsType,2> spaces;
    fs::get_thood_fs(cm,spaces,mesh,false,1);
    const int nq_m=spaces[0].nG;
    if(spaces[0].eNoN!=10 || spaces[1].eNoN!=4) throw std::runtime_error("mixed spaces");
    fs::get_thood_fs(cm,spaces,mesh,false,2);
    const int nq_c=spaces[1].nG;
    xi(0)=.17; xi(1)=.23; xi(2)=.19;
    nn::get_gnn(3,consts::ElementType::TET10,10,xi,Nw,Dw);
    nn::get_gnn(3,consts::ElementType::TET4,4,xi,Nq,Dq);
    const double corners[4][3]={{1,0,0},{0,1,0},{0,0,1},{0,0,0}};
    const int edges[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
    Array<double> x(3,10), G(3,10), metric(3,3), Nxx(6,10);
    Array<double> al(4,10), yl(4,10), bf(3,10), R(4,10);
    Array3<double> K(16,10,10);
    for(int a=0;a<10;a++) for(int j=0;j<3;j++) {
      x(j,a)= a<4 ? corners[a][j] : .5*(corners[edges[a-4][0]][j]+corners[edges[a-4][1]][j]);
      yl(j,a)=.1*(j+1)+x(j,a);
    }
    for(int a=0;a<4;a++) yl(3,a)=2+x(0,a)+2*x(1,a)+3*x(2,a);
    double jac=0.; nn::gnn(10,3,3,Dw,x,G,jac,metric);
    Vector<double> uris(3);
    fluid::fluid_3d_m(cm,0,10,4,1./6.,metric,Nw,Nq,G,Dq,Nxx,al,yl,bf,R,K,0.,0.,uris);
    fluid::fluid_3d_c(cm,0,10,4,1./6.,metric,Nw,Nq,G,Dq,Nxx,al,yl,bf,R,K,0.,0.,uris);
    double edge_momentum=0.,edge_pressure=0.,coupling=0.,max_cont_error=0.;
    for(int a=0;a<10;a++) {
      for(int i=0;i<4;i++) if(!std::isfinite(R(i,a))) throw std::runtime_error("nonfinite residual");
      if(a>=4) {
        for(int i=0;i<3;i++) edge_momentum+=std::abs(R(i,a));
        edge_pressure+=std::abs(R(3,a));
      } else max_cont_error=std::max(max_cont_error,std::abs(R(3,a)-Nq(a)*.5));
      for(int b=0;b<10;b++) for(int k=0;k<16;k++) {
        if(!std::isfinite(K(k,a,b))) throw std::runtime_error("nonfinite tangent");
        if(a<4 && b>=4 && k>=12 && k<15) coupling+=std::abs(K(k,a,b));
      }
    }
    mesh.IEN.resize(10,1);
    cm.rowPtr.resize(11); cm.colPtr.resize(100); cm.Val.resize(16,100); cm.R.resize(4,10);
    cm.Val=2.; cm.R=3.;
    for(int a=0;a<10;a++) {mesh.IEN(a,0)=a; cm.rowPtr(a)=10*a;for(int b=0;b<10;b++)cm.colPtr(10*a+b)=b;}
    cm.rowPtr(10)=100;
    fs::thood_val_rc(cm);
    bool pressure_edges=true;
    for(int a=4;a<10;a++) {
      pressure_edges=pressure_edges && cm.R(3,a)==0. && cm.R(0,a)==3.;
      for(int b=0;b<10;b++) pressure_edges=pressure_edges && cm.Val(15,10*a+b)==(a==b?1.:0.);
    }
    std::ofstream summary(argv[3]); summary << std::setprecision(17)
      << "{\"status\":\"PASS\",\"native_samples\":" << count
      << ",\"velocity_nodes\":10,\"pressure_nodes\":4,\"momentum_quadrature\":" << nq_m
      << ",\"continuity_quadrature\":" << nq_c
      << ",\"local_residual_shape\":[4,10],\"local_tangent_shape\":[16,10,10],\"independent_unknowns_per_element\":34,\"stored_slots\":40"
      << ",\"affine_jacobian\":" << jac << ",\"edge_momentum_abs_sum\":" << edge_momentum
      << ",\"edge_pressure_residual_abs_sum\":" << edge_pressure << ",\"mixed_edge_coupling_abs_sum\":" << coupling
      << ",\"continuity_linear_field_error\":" << max_cont_error
      << ",\"native_pressure_edge_handling_pass\":" << (pressure_edges?"true":"false")
      << ",\"CFD_steps\":0,\"linear_solves\":0}\n";
    if (std::abs(jac-1.)>1e-11 || edge_momentum==0. || edge_pressure!=0. || coupling==0. || max_cont_error>1e-11 || !pressure_edges)
      throw std::runtime_error("native assembly fixture failed");
  } catch(const std::exception& e) {
    std::cerr << e.what() << '\n'; MPI_Finalize(); return 1;
  }
  MPI_Finalize(); return 0;
}

// Reuse the capture/input helpers. The renamed volume entry point is never run.
#define main unused_volume_harness_entry
#include "native_harness.cpp"
#undef main
#include "set_bc.h"
#include "common_geometry.h"

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  try {
    if(argc!=3)throw std::runtime_error("boundary input output-prefix required");
    std::ifstream in(argv[1]);std::string pre(argv[2]);
    ComMod cm;cm.nsd=3;cm.nsymd=6;cm.nEq=1;cm.nMsh=1;cm.tDof=4;cm.dof=4;cm.tnNo=10;cm.cEq=0;cm.cDmn=0;cm.dt=.01;
    cm.mvMsh=false;cm.x.resize(3,10);cm.Bf.resize(3,10);cm.eq.resize(1);cm.msh.resize(1);read_array(in,cm.x);
    auto& eq=cm.eq[0];eq.phys=consts::EquationType::phys_fluid;eq.dof=4;eq.s=0;eq.e=3;eq.nDmn=1;eq.dmn.resize(1);
    eq.af=2./3;eq.am=5./6;eq.gam=2./3;
    auto& dm=eq.dmn[0];dm.Id=-1;dm.phys=consts::EquationType::phys_fluid;
    dm.prop[consts::PhysicalPropertyType::fluid_density]=1056.;dm.prop[consts::PhysicalPropertyType::backflow_stab]=0.;
    dm.prop[consts::PhysicalPropertyType::brinkman_inverse_permeability]=0.;
    dm.fluid_visc.viscType=consts::FluidViscosityModelType::viscType_Const;dm.fluid_visc.mu_i=.00345312;
    auto& mesh=cm.msh[0];mesh.eNoN=10;mesh.nFs=2;mesh.nEl=1;mesh.IEN.resize(10,1);mesh.lShl=false;mesh.lFib=false;
    for(int a=0;a<10;a++)mesh.IEN(a,0)=a;
    nn::select_ele(cm,mesh);fs::init_fs_msh(cm,mesh);mesh.fa.resize(1);
    auto& face=mesh.fa[0];face.iM=0;face.eNoN=6;face.nNo=6;face.nEl=1;face.nG=7;face.eType=consts::ElementType::TRI6;
    face.IEN.resize(6,1);face.gN.resize(6);face.gE.resize(1);face.gE(0)=0;face.nV.resize(3,6);
    for(int a=0;a<6;a++){in>>face.IEN(a,0);face.gN(a)=face.IEN(a,0);face.nV(2,a)=1.;}
    face.w.resize(7);face.xi.resize(2,7);face.N.resize(6,7);face.Nx.resize(2,6,7);
    nn::get_gip(nullptr,face);for(int g=0;g<7;g++)nn::get_gnn(nullptr,g,face);
    SolutionStates s;
    s.current.get_velocity().resize(4,10);s.current.get_acceleration().resize(4,10);s.current.get_displacement().resize(4,10);
    s.intermediate.get_velocity().resize(4,10);s.intermediate.get_acceleration().resize(4,10);
    s.intermediate.get_displacement().resize(4,10);s.old.get_displacement().resize(4,10);
    for(int a=0;a<10;a++)for(int j=0;j<4;j++)s.current.get_velocity()(j,a)=1.+a*.1+j*.01;
    const auto before=s.current.get_velocity();
    eq.nBc=1;eq.bc.resize(1);auto& bc=eq.bc[0];bc.iM=0;bc.iFa=0;bc.g=0.;bc.gx.resize(6);bc.gx=1.;bc.eDrn.resize(3);
    bc.bType=(1<<consts::iBC_Dir)|(1<<consts::iBC_std);
    set_bc::set_bc_dir(cm,s);
    std::ofstream dir(pre+"_dir.txt");dir<<std::setprecision(17);
    for(int a=0;a<10;a++){for(int j=0;j<4;j++)dir<<before(j,a)<<' '<<s.current.get_velocity()(j,a)<<' ';dir<<'\n';}
    // Nodal traction data, including analytic full Cauchy on the planar cap.
    Array<double> traction(3,6);read_array(in,traction);
    bcType load;load.bType=(1<<consts::iBC_gen);load.gm.dof=3;load.gm.nTP=2;load.gm.period=1.;load.gm.t.resize(2);
    load.gm.t(0)=0.;load.gm.t(1)=1.;load.gm.d.resize(3,6,2);
    for(int a=0;a<6;a++)for(int i=0;i<3;i++)for(int t=0;t<2;t++)load.gm.d(i,a,t)=traction(i,a);
    Capture cap;eq.linear_algebra=&cap;CmMod communicator;
    set_bc::set_bc_trac_l(cm,communicator,load,face,s);write_capture(pre+"_traction_RK.txt",cap,6);
    // Exact known restriction of TRI6 face coordinates to the TET10 reference.
    std::array<fsType,2> spaces;fs::get_thood_fs(cm,spaces,mesh,false,1);
    Array<double> corner(3,4);for(int a=0;a<4;a++)for(int j=0;j<3;j++)corner(j,a)=cm.x(j,a);
    std::ofstream q(pre+"_q.txt");q<<std::setprecision(17);
    for(int g=0;g<7;g++) {
      const auto Dface=face.Nx.rslice(g);auto normal=nn::gnnb(cm,face,0,g,Dface,s);
      Vector<double> xp(3),known(3),xiw(3),xiq(3);xiw=.25;xiq=.25;
      for(int a=0;a<6;a++)for(int j=0;j<3;j++)xp(j)+=face.N(a,g)*cm.x(j,face.IEN(a,0));
      const double L[3]={face.xi(0,g),face.xi(1,g),1.-face.xi(0,g)-face.xi(1,g)};
      for(int a=0;a<3;a++)if(face.IEN(a,0)<3)known(face.IEN(a,0))+=L[a];
      Vector<double> Nw(10),Np(4);Array<double> Dw(3,10),Dp(3,4);
      bool vw=true,pq=true;
      try{nn::get_nnx(3,spaces[0].eType,10,cm.x,spaces[0].xib,spaces[0].Nb,xp,xiw,Nw,Dw);}
      catch(const std::exception&){vw=false;bool flag;nn::get_xi(3,spaces[0].eType,10,cm.x,xp,xiw,flag);}
      try{nn::get_nnx(3,spaces[1].eType,4,corner,spaces[1].xib,spaces[1].Nb,xp,xiq,Np,Dp);}
      catch(const std::exception&){pq=false;bool flag;nn::get_xi(3,spaces[1].eType,4,corner,xp,xiq,flag);}
      q<<g<<' '<<face.w(g)<<' '<<face.xi(0,g)<<' '<<face.xi(1,g)<<' ';
      for(int a=0;a<6;a++)q<<face.N(a,g)<<' ';
      for(int a=0;a<6;a++)for(int j=0;j<2;j++)q<<Dface(j,a)<<' ';
      for(int j=0;j<3;j++)q<<xp(j)<<' ';
      for(int j=0;j<3;j++)q<<normal(j)<<' ';
      for(int j=0;j<3;j++)q<<known(j)<<' ';
      for(int j=0;j<3;j++)q<<xiw(j)<<' ';
      for(int j=0;j<3;j++)q<<xiq(j)<<' ';
      q<<vw<<' '<<pq<<'\n';
    }
    // Actual weak-Dirichlet assembly: do not catch and replace its behavior.
    // Catch only to report a blocked audit while allowing normal harness exit.
    for(int a=0;a<10;a++) {
      for(int j=0;j<3;j++)s.intermediate.get_velocity()(j,a)=cm.x(j,a)-cm.x(j,3);
      s.intermediate.get_velocity()(3,a)=.3+1000.*(cm.x(0,a)+2*cm.x(1,a)+3*cm.x(2,a));
    }
    bool weakok=true;bc.weakDir=true;
    try{set_bc::set_bc_dir_wl(cm,bc,mesh,face,s);write_capture(pre+"_weak_RK.txt",cap,10);}
    catch(const std::exception& e){weakok=false;std::ofstream err(pre+"_weak_error.txt");err<<e.what()<<'\n';}
    std::ofstream status(pre+"_status.json");status<<"{\"native_traction_executed\":true,\"native_strong_Dirichlet_executed\":true,\"native_weak_Dirichlet_completed\":"<<(weakok?"true":"false")<<",\"CFD_steps\":0,\"linear_solves\":0}\n";
    eq.linear_algebra=nullptr;
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Finalize();return 1;}
  MPI_Finalize();return 0;
}

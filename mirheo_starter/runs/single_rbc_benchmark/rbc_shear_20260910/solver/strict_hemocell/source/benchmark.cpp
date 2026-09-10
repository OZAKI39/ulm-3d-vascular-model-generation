// Derived from HemoCell oneCellShear/stretchCell and the local pure-fluid adapter.
// Copyright HemoCell contributors, University of Amsterdam. AGPL-3.0-or-later.
// See COPYING. Native library and official examples are not modified.
#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "palabos3D.h"
#include "palabos3D.hh"
#include <mpi.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <memory>
#include <unistd.h>
using namespace hemo;
using Clock=std::chrono::steady_clock;
double sec(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
double maxrank(double x){double y;MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);return y;}

int main(int argc,char** argv){
 if(argc!=2 || std::string(argv[1])=="--help"){
  std::cout<<"Usage: single_rbc_shear_benchmark config.xml (CPU/MPI; fresh case only)\n";return argc==2?0:2;
 }
 auto start=Clock::now(); HemoCell hc(argv[1],argc,argv); auto &cfg=*hc.cfg;
 int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 param::lbm_base_parameters(cfg);
 const double dx=cfg["benchmark"]["dxStar"].read<T>(), dt=cfg["benchmark"]["dtStar"].read<T>();
 const double L=cfg["benchmark"]["length"].read<T>(), H=cfg["benchmark"]["gap"].read<T>();
 const double shear=cfg["benchmark"]["shear"].read<T>(), energy=cfg["domain"]["kBT"].read<T>();
 const int nx=std::lround(L/dx),ny=nx,nz=std::lround(H/dx)+1;
 const int hasCell=cfg["benchmark"]["cell"].read<int>(),probe=cfg["benchmark"]["probe"].read<int>();
 const int prep=cfg["sim"]["prepSteps"].read<int>(),steps=cfg["sim"]["steps"].read<int>(),every=cfg["sim"]["sampleEvery"].read<int>();
 if(nx<8||nz<8||dt<=0||every<1||steps<0||prep<0||hasCell<0||hasCell>1||param::tau<=.5)throw std::runtime_error("INVALID_CASE");
 hc.lattice=new MultiBlockLattice3D<T,DESCRIPTOR>(defaultMultiBlockPolicy3D().getMultiBlockManagement(nx,ny,nz,2),
  defaultMultiBlockPolicy3D().getBlockCommunicator(),defaultMultiBlockPolicy3D().getCombinedStatistics(),
  defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(1./param::tau));
 auto &lat=*hc.lattice;lat.toggleInternalStatistics(false);
 lat.periodicity().toggle(0,true);lat.periodicity().toggle(1,true);lat.periodicity().toggle(2,false);
 std::unique_ptr<OnLatticeBoundaryCondition3D<T,DESCRIPTOR>> bc(createLocalBoundaryCondition3D<T,DESCRIPTOR>());
 Box3D bottom(0,nx-1,0,ny-1,0,0),top(0,nx-1,0,ny-1,nz-1,nz-1);
 bc->setVelocityConditionOnBlockBoundaries(lat,bottom);bc->setVelocityConditionOnBlockBoundaries(lat,top);
 auto walls=[&](double rate){double v=rate*H*.5*dt/dx;
  setBoundaryVelocity(lat,bottom,plb::Array<T,3>(-v,0.,0.));setBoundaryVelocity(lat,top,plb::Array<T,3>(v,0.,0.));};
 initializeAtEquilibrium(lat,lat.getBoundingBox(),1.,plb::Array<T,3>(0.,0.,0.));
 setExternalVector(lat,lat.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,plb::Array<T,3>(0.,0.,0.));walls(0);lat.initialize();
 HemoCellField *cf=nullptr;int nv=0;
 if(hasCell){hc.initializeCellfield();hc.addCellType<RbcHighOrderModel>("RBC",RBC_FROM_SPHERE);
  hc.setMaterialTimeScaleSeparation("RBC",1);hc.setParticleVelocityUpdateTimeScaleSeparation(1);
  hc.loadParticles();cf=(*hc.cellfields)["RBC"];nv=cf->numVertex;
  if(hc.cellfields->number_of_cells!=1)throw std::runtime_error("EXACTLY_ONE_CELL_REQUIRED");
  hc.cellfields->syncEnvelopes();hc.cellfields->applyConstitutiveModel(true);
 }
 auto vertices=[&](){
  std::vector<double> local(nv*7,0),all(nv*7,0);
  if(hasCell){auto &pf=*cf->getParticleField3D();auto &pm=pf.getMultiBlockManagement();
   for(auto id:pm.getLocalInfo().getBlocks()){
    auto b=pm.getUniqueBulk(id);auto &block=pf.getComponent(id);
    for(auto &p:block.particles){auto x=p.sv.position;
     if(!std::isfinite(x[0]+x[1]+x[2]))throw std::runtime_error("NONFINITE_MEMBRANE_COORDINATE");
     if(x[0]<b.x0-.5||x[0]>=b.x1+.5||x[1]<b.y0-.5||x[1]>=b.y1+.5||x[2]<b.z0-.5||x[2]>=b.z1+.5)continue;
     if(hc.cellfields->base_cell_id(p.sv.cellId)!=0||p.sv.vertexId>=nv){
      std::cerr<<"BAD_PARTICLE cell="<<p.sv.cellId<<" vertex="<<p.sv.vertexId<<" expected_vertices="<<nv<<" pos="<<x[0]<<","<<x[1]<<","<<x[2]<<"\n";
      throw std::runtime_error("UNEXPECTED_PARTICLE_ID");}
     int k=p.sv.vertexId*7;local[k]++;
     for(int j=0;j<3;j++){local[k+1+j]+=x[j]*dx;local[k+4+j]+=p.sv.force[j]*param::df/(energy/1e-6);}
    }
   }
   MPI_Allreduce(local.data(),all.data(),all.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   for(int i=0;i<nv;i++)if(all[i*7]!=1){std::cerr<<"VERTEX_OWNERSHIP id="<<i<<" count="<<all[i*7]<<" iter="<<hc.iter<<"\n";
    auto &pf=*cf->getParticleField3D();for(auto id:pf.getMultiBlockManagement().getLocalInfo().getBlocks())for(auto &p:pf.getComponent(id).particles)if(p.sv.vertexId==i)std::cerr<<"CANDIDATE cell="<<p.sv.cellId<<" xyz="<<p.sv.position[0]<<","<<p.sv.position[1]<<","<<p.sv.position[2]<<"\n";
    throw std::runtime_error("MISSING_OR_DUPLICATE_VERTEX");}
  }return all;
 };
 auto init=vertices();
 if(rank==0&&hasCell){std::ofstream mesh("reference.off");mesh<<std::setprecision(17)<<"OFF\n"<<nv<<" "<<cf->triangle_list.size()<<" 0\n";
  for(int i=0;i<nv;i++)mesh<<init[7*i+1]-L/2<<" "<<init[7*i+2]-L/2<<" "<<init[7*i+3]-H/2<<"\n";
  for(auto f:cf->triangle_list)mesh<<"3 "<<f[0]<<" "<<f[1]<<" "<<f[2]<<"\n";
 }
 // Independent affine response: call the actual installed native constitutive model.
 // No surrogate force model, no shear-trajectory fitting, no dynamic claim from this probe.
 if(probe&&hasCell){
  if(rank==0){std::ofstream f("native_response.csv");f<<std::setprecision(17)<<"mode,epsilon,vertex,x,y,z,fx,fy,fz\n";
   for(int mode=0;mode<2;mode++)for(double e:{-.05,-.03,-.01,0.,.01,.03,.05}){
    std::vector<std::unique_ptr<HemoCellParticle>> owned;std::map<int,std::vector<HemoCellParticle*>> cells;std::map<int,bool> present{{0,true}};
    for(int i=0;i<nv;i++){
     double x=init[7*i+1]-L/2,y=init[7*i+2]-L/2,z=init[7*i+3]-H/2;
     if(mode==0)x+=e*z;else{x*=1+e;z/=1+e;}
     owned.emplace_back(new HemoCellParticle({x/dx,y/dx,z/dx},0,i,cf->ctype));cells[0].push_back(owned.back().get());
    }
    cf->mechanics->ParticleMechanics(cells,present,cf->ctype);
    for(int i=0;i<nv;i++){auto &p=*owned[i];f<<(mode==0?"shear":"stretch")<<','<<e<<','<<i;
     for(int j=0;j<3;j++)f<<','<<p.sv.position[j]*dx;
     for(int j=0;j<3;j++)f<<','<<p.sv.force[j]*param::df/(energy/1e-6);f<<'\n';
    }
   }
  }
 }
 std::ofstream membrane,profile,timing,localflow;
 if(rank==0){membrane.open("vertices.csv");profile.open("profiles.csv");timing.open("timings.csv");
  membrane<<std::setprecision(17)<<"step,phase,time_star,strain,vertex,x,y,z,fx,fy,fz\n";
  profile<<std::setprecision(17)<<"step,phase,time_star,strain,bin,z,count,ux,uy,uz,density\n";
  timing<<std::setprecision(17)<<"step,phase,chunk_steps,compute_s,output_s\n";
  localflow.open("local_flow.csv");localflow<<std::setprecision(17)<<"step,phase,time_star,strain,ix,iy,iz,x,y,z,count,ux,uy,uz,density\n";}
 auto sample=[&](int done){
  double t=(done-prep)*dt,g=std::max(0.,t)*shear;std::string phase=done<prep?"relaxation":"shear";auto v=vertices();
  if(rank==0)for(int i=0;i<nv;i++){membrane<<done<<','<<phase<<','<<t<<','<<g<<','<<i;for(int j=1;j<7;j++)membrane<<','<<v[7*i+j];membrane<<'\n';}
  std::vector<double>a(nz*5,0),b(nz*5,0),la(512*5,0),lb(512*5,0);auto &mg=lat.getMultiBlockManagement();
  for(auto id:mg.getLocalInfo().getBlocks()){auto box=mg.getUniqueBulk(id);auto &bl=lat.getComponent(id);auto o=bl.getLocation();
   for(int x=box.x0;x<=box.x1;x++)for(int y=box.y0;y<=box.y1;y++)for(int z=box.z0;z<=box.z1;z++){
    auto &c=bl.get(x-o.x,y-o.y,z-o.z);plb::Array<T,3>u;c.computeVelocity(u);double rho=c.computeDensity();
    if(!std::isfinite(rho+u[0]+u[1]+u[2])||rho<=0)throw std::runtime_error("NONFINITE_FLUID");
    a[5*z]++;for(int j=0;j<3;j++)a[5*z+1+j]+=u[j]*dx/dt;a[5*z+4]+=rho*8;
    // Interior lattice nodes only; local bins share the DPD 8 x 8 x 8 partition.
    if(z>0&&z<nz-1){int ix=std::min(7,int(x*dx/L*8)),iy=std::min(7,int(y*dx/L*8)),iz=std::min(7,int(z*dx/H*8)),k=5*((ix*8+iy)*8+iz);
     la[k]++;for(int j=0;j<3;j++)la[k+1+j]+=u[j]*dx/dt;la[k+4]+=rho*8;}
   }}
  MPI_Allreduce(a.data(),b.data(),b.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(la.data(),lb.data(),lb.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if(rank==0){for(int ix=0;ix<8;ix++)for(int iy=0;iy<8;iy++)for(int iz=0;iz<8;iz++){
   int k=5*((ix*8+iy)*8+iz);localflow<<done<<','<<phase<<','<<t<<','<<g<<','<<ix<<','<<iy<<','<<iz<<','<<(ix+.5)*L/8<<','<<(iy+.5)*L/8<<','<<(iz+.5)*H/8<<','<<lb[k];
   for(int j=1;j<5;j++)localflow<<','<<lb[k+j]/std::max(1.,lb[k]);localflow<<'\n';}localflow.flush();}
  if(rank==0){for(int z=0;z<nz;z++){profile<<done<<','<<phase<<','<<t<<','<<g<<','<<z<<','<<z*dx<<','<<b[5*z];for(int j=1;j<5;j++)profile<<','<<b[5*z+j]/b[5*z];profile<<'\n';}membrane.flush();profile.flush();}
 };
 double setup=maxrank(sec(start)),compute=0,output=0,relax=0;auto s0=Clock::now();sample(0);output+=maxrank(sec(s0));int done=0;
 while(done<prep+steps){
  int stop=rank==0&&access("STOP_REQUESTED",F_OK)==0;MPI_Bcast(&stop,1,MPI_INT,0,MPI_COMM_WORLD);if(stop)break;
  if(done==prep)walls(shear);
  int n=std::min(every,prep+steps-done);if(done<prep)n=std::min(n,prep-done);
  bool relaxing=done<prep;MPI_Barrier(MPI_COMM_WORLD);auto t=Clock::now();
  for(int i=0;i<n;i++){if(hasCell)hc.iterate();else{lat.collideAndStream();hc.iter++;}}
  double c=maxrank(sec(t));if(relaxing)relax+=c;else compute+=c;done+=n;t=Clock::now();sample(done);double o=maxrank(sec(t));output+=o;
  if(rank==0){timing<<done<<','<<(relaxing?"relaxation":"shear")<<','<<n<<','<<c<<','<<o<<'\n';timing.flush();}
 }
 auto elapsed=maxrank(sec(start));
 if(rank==0){std::ofstream o("completion.json");o<<std::setprecision(17)<<"{\"completed\":"<<(done==prep+steps?"true":"false")
  <<",\"cell_count\":"<<hasCell<<",\"vertices\":"<<nv<<",\"actual_steps\":"<<done<<",\"prep_steps\":"<<prep
  <<",\"shear_steps\":"<<std::max(0,done-prep)<<",\"strain_end\":"<<std::max(0,done-prep)*dt*shear<<",\"dt_star\":"<<dt
  <<",\"setup_s\":"<<setup<<",\"relaxation_s\":"<<relax<<",\"coupled_s\":"<<compute<<",\"output_s\":"<<output<<",\"worker_elapsed_s\":"<<elapsed
  <<",\"rank_count\":"<<size<<",\"precision_bits\":"<<sizeof(T)*8<<",\"tau\":"<<param::tau<<",\"lattice_nodes\":"<<nx*ny*nz
  <<",\"material_calls\":"<<(hasCell?done:0)<<",\"IBM_update_interval\":1,\"material_update_interval\":1,\"checkpoints_loaded\":0,\"cellfield_created\":"<<(hasCell?"true":"false")<<"}\n";}
 return 0;
}

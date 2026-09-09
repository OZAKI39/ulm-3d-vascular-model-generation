/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab
in the University of Amsterdam. Any questions or remarks regarding this library
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
// Local derivative: periodic cell-free mean-flow benchmark, 2026.
#include <hemocell.h>
#include "palabos3D.h"
#include "palabos3D.hh"
#include <mpi.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <unistd.h>
using T = double;
using namespace hemo;
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
double maximum(double x) { double v; MPI_Allreduce(&x,&v,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); return v; }
int main(int argc,char** argv) {
 if(argc!=2 || std::string(argv[1])=="--help") { std::cout<<"Usage: pure_fluid_benchmark config.xml\n";return argc==2?0:2; }
 auto start=Clock::now();
 HemoCell hemocell(argv[1],argc,argv);
 int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 auto &cfg=*hemocell.cfg;
 int N=cfg["domain"]["N"].read<int>(), bins=cfg["domain"]["bins"].read<int>();
 double dx=cfg["domain"]["dx"].read<T>(), dt=cfg["domain"]["dt"].read<T>();
 double nu=cfg["domain"]["nuP"].read<T>(),rhoSI=cfg["domain"]["rhoP"].read<T>();
 double acceleration=cfg["domain"]["acceleration"].read<T>();
 int steps=cfg["sim"]["steps"].read<int>(),every=cfg["sim"]["sampleEvery"].read<int>();
 int warm=cfg["sim"]["warmupSteps"].read<int>();
 if(N%2 || N%bins || bins%2 || dx<=0 || dt<=0 || every<=0 || steps<=0 || warm<0 || warm>=steps) throw std::runtime_error("INVALID_CONFIG");
 double tau=.5+3*nu*dt/(dx*dx), a=acceleration*dt*dt/dx;
 param::tau=tau;param::dx=dx;param::dt=dt;
 auto geomStart=Clock::now();
 hemocell.lattice=new MultiBlockLattice3D<T,DESCRIPTOR>(N,N,N,new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(1/tau));
 auto &lat=*hemocell.lattice;
 lat.toggleInternalStatistics(false);lat.periodicity().toggleAll(true);
 double geometryS=maximum(seconds(geomStart));
 initializeAtEquilibrium(lat,lat.getBoundingBox(),1.,plb::Array<T,3>(0.,0.,0.));
 setExternalVector(lat,lat.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,plb::Array<T,3>(0.,0.,0.));
 lat.initialize();
 auto &mgmt=lat.getMultiBlockManagement();
 auto visit=[&](auto action) {
   for(auto id:mgmt.getLocalInfo().getBlocks()) {
    auto b=mgmt.getUniqueBulk(id);auto &block=lat.getComponent(id);auto origin=block.getLocation();
    for(plint x=b.x0;x<=b.x1;++x)for(plint y=b.y0;y<=b.y1;++y)for(plint z=b.z0;z<=b.z1;++z)
     action(block.get(x-origin.x,y-origin.y,z-origin.z),x,y,z);
   }
 };
 auto turnOn=[&]() {
   // At switch populations are rest. Guo observable is j/rho + a/2.
   for(int y=0;y<N;++y) {
    double f=(y<N/2?-a:a);Box3D slab(0,N-1,y,y,0,N-1);
    initializeAtEquilibrium(lat,slab,1.,plb::Array<T,3>(-.5*f,0.,0.));
    setExternalVector(lat,slab,DESCRIPTOR<T>::ExternalField::forceBeginsAt,plb::Array<T,3>(f,0.,0.));
   }
   lat.initialize();
 };
 if(warm==0)turnOn();
 std::ofstream profile,timing;
 if(rank==0) {
  profile.open("profiles.csv"); timing.open("timings.csv");profile<<std::setprecision(17);timing<<std::setprecision(17);
  profile<<"step,time_si,bin,y_si,count,ux_si,uy_si,uz_si,rho_si,rho_min_si,rho_max_si,raw_ux_si,half_force_si,max_speed_si\n";
  timing<<"step,time_si,chunk_steps,compute_max_rank_s,sampling_max_rank_s,compute_cumulative_s,sampling_cumulative_s,elapsed_worker_s\n";
 }
 double computeS=0,sampleS=0,forceSetupS=0;
 auto sample=[&](int step) {
  std::vector<double> local(bins*10,0),global(bins*10,0),lo(bins,1e30),hi(bins,-1e30),maxv(bins,0),glo(bins),ghi(bins),gmv(bins);
  visit([&](Cell<T,DESCRIPTOR>&cell,plint,plint y,plint){
   int j=y*bins/N;plb::Array<T,3>u,momentum;cell.computeVelocity(u);momentTemplates<T,DESCRIPTOR>::get_j(cell,momentum);
   double rho=cell.computeDensity();double f=cell.getExternal(DESCRIPTOR<T>::ExternalField::forceBeginsAt)[0];
   double vnorm=std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
   if(!std::isfinite(rho+vnorm) || rho<=0)MPI_Abort(MPI_COMM_WORLD,3);
   local[j*10]++;for(int k=0;k<3;k++)local[j*10+1+k]+=u[k];local[j*10+4]+=rho;
   local[j*10+5]+=momentum[0]/rho;local[j*10+6]+=f/2;
   lo[j]=std::min(lo[j],rho);hi[j]=std::max(hi[j],rho);maxv[j]=std::max(maxv[j],vnorm);
  });
  MPI_Allreduce(local.data(),global.data(),bins*10,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(lo.data(),glo.data(),bins,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(hi.data(),ghi.data(),bins,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(maxv.data(),gmv.data(),bins,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  double total=0;for(int b=0;b<bins;b++)total+=global[b*10];
  if(total!=N*N*N)MPI_Abort(MPI_COMM_WORLD,4);
  if(rank==0)for(int b=0;b<bins;b++) {
   double n=global[b*10],vscale=dx/dt;
   profile<<step<<','<<step*dt<<','<<b<<','<<(b+.5)*N*dx/bins<<','<<n;
   for(int k=1;k<=3;k++)profile<<','<<global[b*10+k]/n*vscale;
   profile<<','<<global[b*10+4]/n*rhoSI<<','<<glo[b]*rhoSI<<','<<ghi[b]*rhoSI
    <<','<<global[b*10+5]/n*vscale<<','<<global[b*10+6]/n*vscale<<','<<gmv[b]*vscale<<'\n';
  }
  if(rank==0)profile.flush();
 };
 double setupS=maximum(seconds(start));sample(0);
 int done=0;
 while(done<steps) {
  int stop=rank==0 && access("STOP_REQUESTED",F_OK)==0;MPI_Bcast(&stop,1,MPI_INT,0,MPI_COMM_WORLD);if(stop)break;
  if(done==warm && warm>0){auto t=Clock::now();turnOn();forceSetupS+=maximum(seconds(t));}
  int n=std::min(every,steps-done);if(done<warm)n=std::min(n,warm-done);
  MPI_Barrier(MPI_COMM_WORLD);auto t=Clock::now();
  for(int i=0;i<n;i++){lat.collideAndStream();hemocell.iter++;}
  double comp=maximum(seconds(t));computeS+=comp;done+=n;
  t=Clock::now();sample(done);double samp=maximum(seconds(t));sampleS+=samp;
  if(rank==0){timing<<done<<','<<done*dt<<','<<n<<','<<comp<<','<<samp<<','<<computeS<<','<<sampleS<<','<<seconds(start)<<'\n';timing.flush();}
 }
 double elapsed=maximum(seconds(start));
 if(rank==0){
  std::ofstream out("completion.json");out<<std::setprecision(17)
   <<"{\"actual_steps\":"<<done<<",\"actual_time_si\":"<<done*dt<<",\"dt_si\":"<<dt
   <<",\"tau\":"<<tau<<",\"dx_si\":"<<dx<<",\"rho_si\":"<<rhoSI<<",\"nu_si\":"<<nu
   <<",\"rank_count\":"<<size<<",\"lattice_nodes\":"<<N*N*N<<",\"fluid_nodes\":"<<N*N*N
   <<",\"cell_count\":0,\"cellfield_created\":false,\"warmup_steps\":"<<warm
   <<",\"precision_bits\":"<<8*sizeof(T)<<",\"acceleration_lbm\":"<<a
   <<",\"geometry_s\":"<<geometryS<<",\"setup_s\":"<<setupS<<",\"force_setup_s\":"<<forceSetupS
   <<",\"compute_s\":"<<computeS<<",\"sampling_s\":"<<sampleS<<",\"worker_elapsed_s\":"<<elapsed
   <<",\"completed\":"<<(done==steps?"true":"false")<<"}\n";
 }
 return 0;
}

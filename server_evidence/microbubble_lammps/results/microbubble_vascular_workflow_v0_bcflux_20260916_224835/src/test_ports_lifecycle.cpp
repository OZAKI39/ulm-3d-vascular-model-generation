#include "workflow_geometry.hpp"
#include "lifecycle.hpp"
#include "atom.h"
#include "library.h"
#include <mpi.h>
#include <iostream>
#include <cassert>
#include <set>
using namespace rigid;
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rank,nr;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nr);try{
 assert(argc==2);std::string root=argv[1];workflow::Geometry geom(root);frozen::FlowField field(root+"/fields/FROZEN_FLOW_FIELD_V0.h5");
 char s0[]="port_test",s1[]="-screen",s2[]="none",s3[]="-log",s4[]="none";char*args[]={s0,s1,s2,s3,s4};auto*l=new LAMMPS_NS::LAMMPS(5,args,MPI_COMM_WORLD);
 auto cmd=[&](std::string s){l->input->one(s);};cmd("units si");cmd("atom_style sphere");cmd("boundary f f f");cmd("atom_modify map array");cmd("processors 1 1 "+std::to_string(nr));
 cmd("region box block 0.00005 0.00021 0.00003 0.00018 0.00009 0.00019 units box");cmd("create_box 1 box");cmd("pair_style zero 8e-6 full");cmd("pair_coeff * *");cmd("neighbor 4e-7 bin");cmd("neigh_modify every 1 delay 0 check no");cmd("run 0 post no");
 int tests=0;std::set<int>owners;
 for(int k=0;k<4;k++){
  auto&p=geom.ports[k];Vec a=sub(p.center,mul(p.normal,2*field.dx)),b=add(p.center,mul(p.normal,2*field.dx));
  auto crossing=geom.crossing(a,b);assert(crossing.port==k);assert(geom.crossing(b,a).port<0);assert(geom.crossing(a,a).port<0);
  // Directly exercise known cap-crossing segments. These are isolated geometry
  // fixtures, not flow trajectories, not timed inlet injections, not evidence
  // that a simulated injected bubble reached an outlet.
  double radius=.375e-6;assert(geom.wall.segment_safe(a,b,radius+1e-10));int id=101+7*k,type=1;
  assert(lammps_create_atoms(l,1,&id,&type,a.data(),nullptr,nullptr,0)==1);
  auto*at=l->atom;for(int i=0;i<at->nlocal;i++){at->radius[i]=radius;at->rmass[i]=1e-15;for(int d=0;d<3;d++)at->omega[i][d]=0;}
  cmd("run 0 post no");int count=at->nlocal,total=0;MPI_Allreduce(&count,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);assert(total==1);int localowner=count?rank:-1,owner=-1;MPI_Allreduce(&localowner,&owner,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);owners.insert(owner);
  for(int i=0;i<at->nlocal;i++){assert(at->tag[i]==id);for(int d=0;d<3;d++)at->x[i][d]=b[d];}
  cmd("run 0 post no");workflow::remove_atoms(l,{id});assert(at->natoms==0);tests+=4;
 }
 // Side-wall crossing with both endpoints separated from the surface must be
 // rejected by continuous segment safety; a point wall distance alone is not enough.
 Vec c=geom.ports[0].center;c=sub(c,mul(geom.ports[0].normal,2.824899180810653e-6));Vec x=c; x[0]+=5e-6;
 assert(geom.wall.distance(c)>.375e-6);assert(!geom.wall.segment_safe(c,x,.375e-6));tests+=2;
 // Remaining tags survive removal; new tags need not be contiguous.
 int ids[2]={501,509},types[2]={1,1};double xx[6]={c[0],c[1],c[2],c[0],c[1],c[2]-4e-6};assert(lammps_create_atoms(l,2,ids,types,xx,nullptr,nullptr,0)==2);
 for(int i=0;i<l->atom->nlocal;i++){l->atom->radius[i]=.375e-6;l->atom->rmass[i]=1e-15;}cmd("run 0 post no");workflow::remove_atoms(l,{501});assert(l->atom->natoms==1);for(int i=0;i<l->atom->nlocal;i++)assert(l->atom->tag[i]==509);workflow::remove_atoms(l,{509});assert(l->atom->natoms==0);tests+=2;
 // Isolated storage migration fixture: retain tag/radius/omega across processor
 // ownership changes. Its displacement is prescribed, never a flow trajectory.
 int mid=777,type=1;Vec destination=sub(geom.ports[1].center,mul(geom.ports[1].normal,2*field.dx));
 assert(lammps_create_atoms(l,1,&mid,&type,c.data(),nullptr,nullptr,0)==1);
 for(int i=0;i<l->atom->nlocal;i++){l->atom->radius[i]=.375e-6;l->atom->rmass[i]=1e-15;for(int k=0;k<3;k++)l->atom->omega[i][k]=k+1;}
 cmd("run 0 post no");int local_owner=l->atom->nlocal?rank:-1,owner_before=-1,owner_after=-1;MPI_Allreduce(&local_owner,&owner_before,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 for(int i=0;i<l->atom->nlocal;i++)for(int k=0;k<3;k++)l->atom->x[i][k]=destination[k];
 cmd("run 0 post no");local_owner=l->atom->nlocal?rank:-1;MPI_Allreduce(&local_owner,&owner_after,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 assert(nr==1||owner_before!=owner_after);assert(l->atom->natoms==1);
 for(int i=0;i<l->atom->nlocal;i++){assert(l->atom->tag[i]==777&&l->atom->radius[i]==.375e-6);for(int k=0;k<3;k++)assert(l->atom->x[i][k]==destination[k]&&l->atom->omega[i][k]==k+1);}
 workflow::remove_atoms(l,{777});assert(l->atom->natoms==0);tests+=3;
 if(!rank)std::cout<<"{\"status\":\"PASS\",\"tests\":"<<tests<<",\"mpi_ranks\":"<<nr<<",\"port_owners_exercised\":"<<owners.size()<<",\"migration_owner_before\":"<<owner_before<<",\"migration_owner_after\":"<<owner_after<<",\"ports\":[\"INLET_BACKFLOW\",\"OUTLET_0\",\"OUTLET_1\",\"OUTLET_2\"],\"scope\":\"GEOMETRIC_CAP_CROSSING_AND_NATIVE_LAMMPS_LIFECYCLE_FIXTURE\",\"actual_flow_outlet_event\":false,\"disclaimer\":\"NOT EXPERIMENTAL CONCENTRATION\"}\n";
 delete l;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();MPI_Finalize();return 0;
 }catch(std::exception const&e){std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}}

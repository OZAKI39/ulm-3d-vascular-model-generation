#include "lammps.h"
#include "input.h"
#include "modify.h"
#include "fix_passive_flow.hpp"
#include <mpi.h>
#include <exception>
#include <iostream>
using namespace LAMMPS_NS;
static Fix*factory(LAMMPS*l,int n,char**a){return new FixPassiveFlow(l,n,a);}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{auto*l=new LAMMPS(argc,argv,MPI_COMM_WORLD);(*l->modify->fix_map)["sonovue/passive/flow"]=&factory;l->input->file();delete l;}catch(std::exception const&e){std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);}MPI_Finalize();}

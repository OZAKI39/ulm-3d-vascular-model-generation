#include "lammps.h"
#include "library.h"
#include "input.h"
#include "modify.h"
#include "fix_passive_flow.hpp"
#include <mpi.h>
#include <exception>
#include <iostream>
using namespace LAMMPS_NS;
static Fix*factory(LAMMPS*l,int n,char**a){return new FixPassiveFlow(l,n,a);}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{auto*l=new LAMMPS(argc,argv,MPI_COMM_WORLD);Modify::fix_styles().set_plugin("sonovue/passive/flow", &factory);l->input->file();delete l;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();}catch(std::exception const&e){std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);}MPI_Finalize();}

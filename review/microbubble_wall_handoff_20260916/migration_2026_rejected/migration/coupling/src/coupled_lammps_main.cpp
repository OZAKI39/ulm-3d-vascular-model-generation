// External embedding adapter. The upstream LAMMPS main/core and libraries remain unchanged.
#include "lammps.h"
#include "library.h"
#include "input.h"
#include "modify.h"
#include "exceptions.h"
#include "fix_frozen_flow_drag.hpp"
#include <mpi.h>
#include <iostream>
using namespace LAMMPS_NS;
static Fix* make_drag(LAMMPS*lmp,int n,char**a){return new FixFrozenFlowDrag(lmp,n,a);}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{auto*lmp=new LAMMPS(argc,argv,MPI_COMM_WORLD);Modify::fix_styles().set_plugin("frozen/flow/drag", &make_drag);lmp->input->file();delete lmp;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();MPI_Finalize();return 0;}catch(std::exception const&e){std::cerr<<"Coupling adapter error: "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);return 1;}}

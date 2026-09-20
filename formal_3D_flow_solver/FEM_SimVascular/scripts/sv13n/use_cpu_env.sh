#!/bin/sh
export PATH=/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/bin:/usr/bin:/bin
export OMPI_CC=/usr/bin/gcc-12 OMPI_CXX=/usr/bin/g++-12 OMPI_FC=/usr/bin/gfortran-12
export LD_LIBRARY_PATH=/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib:${LD_LIBRARY_PATH:-}
exec "$@"

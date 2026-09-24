#!/usr/bin/env bash
set -euo pipefail
export CUDA_HOME="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3j/external/compat_cuda/cuda-12.3.2"
export PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3j/external/compat_cuda/cuda-12.3.2/bin:/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3j/external/compat_cuda/host-bin:/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3j/external/compat_cuda/cuda-12.3.2/lib64:/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CC=/usr/bin/gcc-12
export CXX=/usr/bin/g++-12
export FC=/usr/bin/gfortran-12
if [[ "$#" -gt 0 ]]; then exec "$@"; fi

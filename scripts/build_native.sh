#!/bin/bash
# Build only. Does not run a CFD or modify the mouse solver installation.
set -euo pipefail
HUMAN_PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SOLVER_DEPS_ROOT=/home/lzy/projects/new_vessel_flow_solvers
BRAVA_NATIVE_ROOT="$SOLVER_DEPS_ROOT/brava_svmultiphysics"
mkdir -p "$BRAVA_NATIVE_ROOT/logs" "$BRAVA_NATIVE_ROOT/bin"
export PATH=/usr/bin:/bin:/usr/local/cuda-12.6/bin
unset PYTHONPATH PYTHONHOME LD_LIBRARY_PATH PETSC_DIR PETSC_ARCH
export OMPI_CC=/usr/bin/gcc-13 OMPI_CXX=/usr/bin/g++-13 OMPI_FC=/usr/bin/gfortran-13
cmake -S "$HUMAN_PROJECT_ROOT/native/svmultiphysics/Code" -B "$BRAVA_NATIVE_ROOT/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc-13 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-13 \
  -DMPI_C_COMPILER=/usr/bin/mpicc -DMPI_CXX_COMPILER=/usr/bin/mpicxx -DMPI_Fortran_COMPILER=/usr/bin/mpifort \
  -DSV_USE_MPI=ON -DSV_USE_DUMMY_MPI=OFF -DSV_USE_TRILINOS=OFF -DENABLE_UNIT_TEST=OFF \
  -DSV_PETSC_DIR="$SOLVER_DEPS_ROOT/petsc-3.25.5-cuda12.6" \
  -DSV_LOCAL_VTK_PATH="$SOLVER_DEPS_ROOT/vtk-9.3/lib/cmake/vtk-9.3" \
  -DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/libblas.so.3 \
  -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/liblapack.so.3 \
  -DCMAKE_INSTALL_PREFIX="$BRAVA_NATIVE_ROOT" > "$BRAVA_NATIVE_ROOT/logs/configure.log" 2>&1
cmake --build "$BRAVA_NATIVE_ROOT/build" --parallel 2 > "$BRAVA_NATIVE_ROOT/logs/build.log" 2>&1
install -m755 "$BRAVA_NATIVE_ROOT/build/bin/svmultiphysics" "$BRAVA_NATIVE_ROOT/bin/svmultiphysics"
echo 'BraVa production source built; no CFD executed.'

#!/usr/bin/env bash
# Source this file to select only the dedicated BraVa CFD environment.
export BRAVA_ROOT=/workspace/brava_cfd
export BRAVA_PROJECT="$BRAVA_ROOT/project"
export BRAVA_SOLVER="$BRAVA_ROOT/build/svmp-cuda12/svMultiPhysics-build/bin/svmultiphysics"
export BRAVA_MPI="$BRAVA_ROOT/bin/brava-mpi"
export VIRTUAL_ENV="$BRAVA_ROOT/envs/cfd"
export PATH="$VIRTUAL_ENV/bin:$BRAVA_ROOT/bin:$BRAVA_ROOT/deps/openmpi-4.1.6/bin:/usr/local/cuda/bin:/usr/local/bin:/usr/bin:/bin"
export OPAL_PREFIX="$BRAVA_ROOT/deps/openmpi-4.1.6"
export OMPI_CC=/usr/bin/gcc OMPI_CXX=/usr/bin/g++
export PETSC_DIR="$BRAVA_ROOT/deps/petsc-cuda12"
export PETSC_ARCH=
export CUDA_HOME=/usr/local/cuda
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$OPAL_PREFIX/lib:$BRAVA_ROOT/deps/vtk-9.3/lib:$BRAVA_ROOT/deps/sysroot/usr/lib/x86_64-linux-gnu:$CUDA_HOME/lib64:/usr/local/nvidia/lib:/usr/local/nvidia/lib64"
export PYTHONPATH="$BRAVA_PROJECT/vendor:$BRAVA_PROJECT/population_v1/source/particle_3d/src"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export UV_CACHE_DIR="$BRAVA_ROOT/runtime/uv-cache"
export CUPY_CACHE_DIR="$BRAVA_ROOT/runtime/cupy-cache"
export CUPY_COMPILE_WITH_PTX=0
export MPLBACKEND=Agg

#!/bin/bash
export Y_CLOT_HOME=/workspace/y_clot
export PATH="$Y_CLOT_HOME/envs/clot/bin:$Y_CLOT_HOME/bin:$Y_CLOT_HOME/deps/openmpi/bin:/usr/local/cuda/bin:$PATH"
export OPAL_PREFIX="$Y_CLOT_HOME/deps/openmpi"
export OMPI_CC=/usr/bin/gcc OMPI_CXX=/usr/bin/g++
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$Y_CLOT_HOME/project/tools/lammps-contact-onset-v1/python:$Y_CLOT_HOME/project:$Y_CLOT_HOME/deps/y_workflow/src:$Y_CLOT_HOME/deps/y_cfd_vendor"
export LD_LIBRARY_PATH="$Y_CLOT_HOME/project/tools/lammps-native-y-wall-shared-v1-build:$Y_CLOT_HOME/deps/petsc/lib:$Y_CLOT_HOME/deps/openmpi/lib:$Y_CLOT_HOME/deps/solver-inputs/vtk/lib:/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export MPI4PY_MPIABI=openmpi MPI4PY_LIBMPI="$Y_CLOT_HOME/deps/openmpi/lib/libmpi.so.40"
export PETSC_DIR="$Y_CLOT_HOME/deps/petsc"
export PYVISTA_OFF_SCREEN=true PYTHONDONTWRITEBYTECODE=1
export UV_CACHE_DIR="$Y_CLOT_HOME/.cache/uv" UV_PYTHON_INSTALL_DIR="$Y_CLOT_HOME/envs/python" UV_PYTHON_BIN_DIR="$Y_CLOT_HOME/bin"

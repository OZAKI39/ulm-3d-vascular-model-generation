#!/usr/bin/env bash

export MIR_ROOT="/home/lzy/projects/mirheo_starter"
export MIR_SRC="$MIR_ROOT/vendor/Mirheo"
export CUDA_HOME="/usr/local/cuda-12.6"

source "$MIR_ROOT/.venv/bin/activate"

export PATH="$VIRTUAL_ENV/bin:$CUDA_HOME/bin:/usr/lib/wsl/lib:/usr/bin:/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:/usr/lib/wsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

export OMPI_CC=/usr/bin/gcc-12
export OMPI_CXX=/usr/bin/g++-12
export OMP_NUM_THREADS=1
export PYTHONNOUSERSITE=1

unset PYTHONPATH
unset PYTHONHOME
unset HDF5_DIR

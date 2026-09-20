#!/usr/bin/env bash
set -euo pipefail
# Stage SV1.3H: scoped CUDA12 environment; no permanent shell or system edits.
export CUDA_HOME="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3h/external/cuda-12.6"
export PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3h/external/cuda-12.6/bin:/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3g/external/gpu_mpi/bin:/usr/bin:/bin${PATH:+:$PATH}"
export LD_LIBRARY_PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3h/external/cuda-12.6/lib64:/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3g/external/gpu_mpi/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
if [[ "$#" -gt 0 ]]; then exec "$@"; fi

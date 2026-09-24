#!/usr/bin/env bash
set -euo pipefail
# Validated Open MPI 4.1.6; no scientific parameters.
export PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3g/external/gpu_mpi/bin:/usr/bin:/bin:/usr/local/cuda/bin"
export LD_LIBRARY_PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3g/external/gpu_mpi/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
root_args=()
if [[ "$(id -u)" == 0 ]]; then root_args=(--allow-run-as-root); fi
exec "/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3g/external/gpu_mpi/bin/mpiexec" "${root_args[@]}" "$@"

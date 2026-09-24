#!/usr/bin/env bash
set -euo pipefail
# Isolated MPI4.1.6 with Fortran bindings; no scientific options.
export PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/bin:${PATH:-/usr/bin:/bin}"
export LD_LIBRARY_PATH="/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
root_args=()
if [[ "$(id -u)" == 0 ]]; then root_args=(--allow-run-as-root); fi
exec "/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/bin/mpiexec" "${root_args[@]}" "$@"

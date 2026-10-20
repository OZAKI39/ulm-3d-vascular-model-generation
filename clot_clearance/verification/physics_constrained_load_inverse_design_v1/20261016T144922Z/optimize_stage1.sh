#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /workspace/y_clot/envs/inverse_opt_v1/bin/python /workspace/y_clot/project/analysis/physics_constrained_load_inverse_design_v1/run_inverse_design.py --limit 2

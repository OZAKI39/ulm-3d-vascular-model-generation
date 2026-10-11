#!/bin/bash
set -e
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec python -m analysis.traction_failure_v1.offline_a

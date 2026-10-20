#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project
exec python -u analysis/clot_shear_load_inverse_design_v1/postprocess.py

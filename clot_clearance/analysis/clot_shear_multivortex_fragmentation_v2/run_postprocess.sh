#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project
exec python -u /workspace/y_clot/project/analysis/clot_shear_multivortex_fragmentation_v2/postprocess.py

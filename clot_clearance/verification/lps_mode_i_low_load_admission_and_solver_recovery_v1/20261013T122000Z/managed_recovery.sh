#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec python -u analysis/mode_i_low_load_v1/recover_saved.py

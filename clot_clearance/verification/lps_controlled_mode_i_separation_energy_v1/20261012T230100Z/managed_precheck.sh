#!/bin/bash
set -e
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
exec timeout --signal=TERM --kill-after=15 340 python -u analysis/mode_i_v1/precheck.py

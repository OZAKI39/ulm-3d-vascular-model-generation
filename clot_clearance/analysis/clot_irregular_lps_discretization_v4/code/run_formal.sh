#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_irregular_lps_discretization_v4/20261019T221000Z/build:${LD_LIBRARY_PATH:-}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project/analysis/clot_irregular_lps_discretization_v4/code
python -u run_case.py pilot
python -u run_case.py continue

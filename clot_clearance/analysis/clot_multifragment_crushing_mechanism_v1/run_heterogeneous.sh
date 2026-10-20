#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_multifragment_crushing_mechanism_v1/20261016T222934Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
exec python -u /workspace/y_clot/project/analysis/clot_multifragment_crushing_mechanism_v1/run_case.py --name heterogeneous --heterogeneous

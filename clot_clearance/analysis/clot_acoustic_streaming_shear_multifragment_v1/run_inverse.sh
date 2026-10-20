#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_multifragment_breakup_v1/20261017T110849Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
exec python -u /workspace/y_clot/project/analysis/clot_acoustic_streaming_shear_multifragment_v1/run_shear_cases.py --name inverse

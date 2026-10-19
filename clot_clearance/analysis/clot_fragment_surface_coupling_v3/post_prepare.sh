#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
exec python -u /workspace/y_clot/project/analysis/clot_fragment_surface_coupling_v3/wait_and_postprocess.py prepare

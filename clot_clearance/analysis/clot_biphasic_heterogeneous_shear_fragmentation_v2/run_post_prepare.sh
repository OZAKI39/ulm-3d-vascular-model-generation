#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYVISTA_OFF_SCREEN=true
cd /workspace/y_clot/project/analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2
exec python -u postprocess.py prepare

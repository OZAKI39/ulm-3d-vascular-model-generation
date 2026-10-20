#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project
exec python -u analysis/clot_fragment_surface_coupling_v3_red_visuals/render.py

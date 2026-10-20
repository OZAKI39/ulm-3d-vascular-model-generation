#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project
python analysis/y_clot_v4_reference_style_v1/code/render.py results
python analysis/y_clot_v4_reference_style_v1/code/render.py orbit

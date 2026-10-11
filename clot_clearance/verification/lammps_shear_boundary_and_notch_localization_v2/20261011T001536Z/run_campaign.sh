#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project
exec python -u /workspace/y_clot/project/verification/lammps_shear_boundary_and_notch_localization_v2/20261011T001536Z/run_campaign.py

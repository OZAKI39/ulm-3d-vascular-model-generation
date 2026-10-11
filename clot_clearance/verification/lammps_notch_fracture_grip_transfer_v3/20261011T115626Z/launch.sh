#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec python -u /workspace/y_clot/project/verification/lammps_notch_fracture_grip_transfer_v3/20261011T115626Z/campaign.py

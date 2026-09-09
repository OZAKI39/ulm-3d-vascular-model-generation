#!/usr/bin/env bash
set -euo pipefail
source /home/lzy/projects/mirheo_starter/scripts/activate_mirheo.sh
exec /usr/bin/mpirun.openmpi --bind-to none -np 2 /home/lzy/projects/mirheo_starter/.venv/bin/python -B -u /home/lzy/projects/mirheo_starter/runs/fluid_calibration/dpd_round1_20260908/high_gamma_equilibrium/gpu_worker.py --spec /home/lzy/projects/mirheo_starter/runs/fluid_calibration/dpd_round1_20260908/high_gamma_equilibrium/actual_parameters.json

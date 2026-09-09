#!/usr/bin/env bash
set -euo pipefail
source /home/lzy/projects/mirheo_starter/scripts/activate_mirheo.sh
exec /usr/bin/mpirun.openmpi --bind-to none -np 2 /home/lzy/projects/mirheo_starter/.venv/bin/python -B -u /home/lzy/projects/mirheo_starter/runs/sdpd_diagnostics/thermal_cause_20260909/equilibrium_half_dt_trend/gpu_worker.py --spec /home/lzy/projects/mirheo_starter/runs/sdpd_diagnostics/thermal_cause_20260909/equilibrium_half_dt_trend/actual_parameters.json

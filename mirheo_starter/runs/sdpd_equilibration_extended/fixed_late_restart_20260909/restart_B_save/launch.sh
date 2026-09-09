#!/usr/bin/env bash
set -euo pipefail
source /home/lzy/projects/mirheo_starter/scripts/activate_mirheo.sh
export PYTHONPATH=/home/lzy/projects/mirheo_starter/runs/sdpd_equilibration_extended/fixed_late_restart_20260909/restart_B_save/frozen_code
exec /usr/bin/mpirun.openmpi --bind-to none -np 2 /home/lzy/projects/mirheo_starter/.venv/bin/python -B -u -m py_scripts.sdpd_diagnostics.extended_worker --spec /home/lzy/projects/mirheo_starter/runs/sdpd_equilibration_extended/fixed_late_restart_20260909/restart_B_save/actual_parameters.json

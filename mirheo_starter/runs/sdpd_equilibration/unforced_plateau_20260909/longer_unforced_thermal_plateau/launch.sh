#!/usr/bin/env bash
set -euo pipefail
source /home/lzy/projects/mirheo_starter/scripts/activate_mirheo.sh
export PYTHONPATH=/home/lzy/projects/mirheo_starter/runs/sdpd_equilibration/unforced_plateau_20260909/longer_unforced_thermal_plateau/frozen_code
exec /usr/bin/mpirun.openmpi --bind-to none -np 2 /home/lzy/projects/mirheo_starter/.venv/bin/python -B -u -m py_scripts.sdpd_diagnostics.equilibration_worker --spec /home/lzy/projects/mirheo_starter/runs/sdpd_equilibration/unforced_plateau_20260909/longer_unforced_thermal_plateau/actual_parameters.json

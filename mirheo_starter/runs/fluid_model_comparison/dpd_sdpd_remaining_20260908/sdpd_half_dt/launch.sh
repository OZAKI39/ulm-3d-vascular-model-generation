#!/usr/bin/env bash
set -euo pipefail
source /home/lzy/projects/mirheo_starter/scripts/activate_mirheo.sh
exec /usr/bin/mpirun.openmpi --bind-to none -np 2 /home/lzy/projects/mirheo_starter/.venv/bin/python -B -u /home/lzy/projects/mirheo_starter/runs/fluid_model_comparison/dpd_sdpd_remaining_20260908/sdpd_half_dt/gpu_worker.py --spec /home/lzy/projects/mirheo_starter/runs/fluid_model_comparison/dpd_sdpd_remaining_20260908/sdpd_half_dt/actual_parameters.json

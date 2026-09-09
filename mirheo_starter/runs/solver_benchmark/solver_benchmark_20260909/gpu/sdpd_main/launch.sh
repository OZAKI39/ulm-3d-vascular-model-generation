#!/bin/bash
set -euo pipefail
source /home/lzy/projects/mirheo_starter/scripts/activate_mirheo.sh
export SOLVER_BENCHMARK_AUTHORIZED=93b1f6561eb339f087476bb17498aba0cc1325d36244e05735eb70f610e09e56
exec /usr/bin/mpirun.openmpi --bind-to none -np 2 /home/lzy/projects/mirheo_starter/.venv/bin/python -B -u /home/lzy/projects/mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main/mirheo_worker.py --spec /home/lzy/projects/mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main/task.json

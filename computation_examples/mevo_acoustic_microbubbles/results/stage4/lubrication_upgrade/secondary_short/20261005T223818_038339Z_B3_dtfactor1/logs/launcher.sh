#!/bin/bash
set -euo pipefail
cd /workspace/mevo_acoustic_microbubbles_stage4l
export PYTHONPATH=src:.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /workspace/mevo_acoustic_microbubbles/.venv/bin/python -u - <<'RUN'
from pathlib import Path
from scripts.stage4l.branch_runs import run_case
print(run_case(Path('/workspace/mevo_acoustic_microbubbles'), 'B3', secondary=True, duration=0.1, max_wall_s=3600), flush=True)
RUN

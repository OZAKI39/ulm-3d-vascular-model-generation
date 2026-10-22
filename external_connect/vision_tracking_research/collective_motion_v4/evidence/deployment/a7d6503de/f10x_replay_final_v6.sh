#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg
export F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10X_SERVER=/workspace/f10x_piv_20261017
export F10X_OUTPUT=/workspace/f10x_piv_20261017/runs/replay_final_v6
source /workspace/f10x_piv_20261017/.venv/bin/activate
cd /workspace/f10x_piv_20261017/snapshots/a7d6503de/external_connect/vision_tracking_research/collective_motion_v4
while supervisorctl status f10x_verified_pipeline_v5b | rg -q RUNNING; do sleep 1; done; python -m pytest -q tests && python scripts/causal_replay.py --strategy regular && python scripts/causal_replay.py --strategy sparse3 && python scripts/causal_replay.py --strategy sparse6

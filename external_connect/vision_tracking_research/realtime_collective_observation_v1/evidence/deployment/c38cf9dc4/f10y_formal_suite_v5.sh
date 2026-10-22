#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg
export F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10X_SERVER=/workspace/f10x_piv_20261017
export F10X_OUTPUT=/workspace/f10y_vision_20261017/runs/formal_v5
export PATH=/workspace/f10y_vision_20261017/.venv/bin:$PATH
cd /workspace/f10y_vision_20261017/snapshots/c38cf9dc4/external_connect/vision_tracking_research/realtime_collective_observation_v1
bash scripts/formal_suite.sh

#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg
export F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10X_SERVER=/workspace/f10x_piv_20261017
export F10X_OUTPUT=/workspace/f10y_vision_20261017/runs/stability_v3
export PATH=/workspace/f10y_vision_20261017/.venv/bin:$PATH
cd /workspace/f10y_vision_20261017/snapshots/e91b1caab/external_connect/vision_tracking_research/realtime_collective_observation_v1
QT_QPA_PLATFORM=offscreen python run_replay.py --backend B4 --video /workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete --autostart --loop --seconds 605 --output /workspace/f10y_vision_20261017/runs/stability_v3

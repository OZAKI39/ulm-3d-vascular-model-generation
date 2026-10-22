#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg
export F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10X_SERVER=/workspace/f10x_piv_20261017
export F10X_OUTPUT=/workspace/f10y_vision_20261017/runs/ui_final_v8
export PATH=/workspace/f10y_vision_20261017/.venv/bin:$PATH
cd /workspace/f10y_vision_20261017/snapshots/2b3707bfa/external_connect/vision_tracking_research/realtime_collective_observation_v1
while supervisorctl status f10y_formal_suite_v6 | grep -q RUNNING; do sleep 5; done
QT_QPA_PLATFORM=offscreen python run_replay.py --backend B4 --video /workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete --autostart --seconds 33 --simulate-events --output /workspace/f10y_vision_20261017/runs/ui_final_v8 --capture-ui /workspace/f10y_vision_20261017/runs/Video_F10Y_1_realtime_replay_observation_v8.mp4

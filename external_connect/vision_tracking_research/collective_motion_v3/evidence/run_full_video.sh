#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export F10W_OUTPUT_ROOT=/workspace/f10w_vision_20261017/runs/full_video_e15b
export F10W_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
mkdir -p "$F10W_OUTPUT_ROOT/results" "$F10W_OUTPUT_ROOT/evidence"
cd /workspace/f10w_vision_20261017/snapshots/e15b0d230/external_connect/vision_tracking_research/collective_motion_v3
/workspace/f10w_vision_20261017/.venv/bin/python scripts/audit_source.py
exec /workspace/f10w_vision_20261017/.venv/bin/python scripts/full_video.py

#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 F10V_DATA_ROOT=/workspace/f10v_vision_20261017/data F10V_OUTPUT_ROOT=/workspace/f10v_vision_20261017/runs/prototype_7aed F10V_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
cd /workspace/f10v_vision_20261017/snapshots/a5a4ef53d/external_connect/vision_tracking_research/collective_motion_v2
/workspace/f10v_vision_20261017/.venv/bin/python scripts/synthetic_demo.py

#!/bin/bash
set -euo pipefail
source /workspace/f10w_vision_20261017/.venv/bin/activate
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
export F10W_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10W_OUTPUT_ROOT=/workspace/f10w_vision_20261017/runs/models_v1
export F10W_THIRD_PARTY=/workspace/f10w_vision_20261017/third_party
export F10W_MODELS=/workspace/f10w_vision_20261017/models
cd /workspace/f10w_vision_20261017/snapshots/7cc85d4f2/external_connect/vision_tracking_research/collective_motion_v3
python scripts/causal_replay.py
python scripts/extra_diagnostics.py real
python scripts/extra_diagnostics.py historical
python scripts/extra_diagnostics.py ot
python -m pip install flowiz
python scripts/extra_diagnostics.py official

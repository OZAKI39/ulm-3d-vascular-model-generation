#!/bin/bash
set -euo pipefail
source /workspace/f10w_vision_20261017/.venv/bin/activate
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
export F10W_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10W_OUTPUT_ROOT=/workspace/f10w_vision_20261017/runs/models_v1
export F10W_THIRD_PARTY=/workspace/f10w_vision_20261017/third_party
export F10W_MODELS=/workspace/f10w_vision_20261017/models
cd /workspace/f10w_vision_20261017/snapshots/2251c5b49/external_connect/vision_tracking_research/collective_motion_v3
mkdir -p "$F10W_OUTPUT_ROOT/evidence/unlite_compat"
cp "$F10W_OUTPUT_ROOT/evidence/model_smoke_UnLiteFlowNet.json" "$F10W_OUTPUT_ROOT/evidence/unlite_compat/attempt1.json"
python scripts/compat_unlite.py
python scripts/model_smoke.py UnLiteFlowNet

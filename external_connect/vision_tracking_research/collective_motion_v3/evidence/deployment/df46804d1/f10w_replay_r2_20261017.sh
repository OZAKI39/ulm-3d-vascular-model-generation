#!/bin/bash
set -euo pipefail
source /workspace/f10w_vision_20261017/.venv/bin/activate
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
export F10W_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10W_OUTPUT_ROOT=/workspace/f10w_vision_20261017/runs/models_v1
export F10W_THIRD_PARTY=/workspace/f10w_vision_20261017/third_party
export F10W_MODELS=/workspace/f10w_vision_20261017/models
cd /workspace/f10w_vision_20261017/snapshots/df46804d1/external_connect/vision_tracking_research/collective_motion_v3
while supervisorctl status f10w_replay_diagnostics_20261017 | grep -q RUNNING; do sleep 2; done
python -m unittest discover -s tests -v
mkdir -p "$F10W_OUTPUT_ROOT/results/replay_before_hash_cache_fix"
for p in "$F10W_OUTPUT_ROOT"/results/causal_replay_*.json; do mv "$p" "$F10W_OUTPUT_ROOT/results/replay_before_hash_cache_fix/"; done
python scripts/causal_replay.py
python scripts/remote_asset_audit.py

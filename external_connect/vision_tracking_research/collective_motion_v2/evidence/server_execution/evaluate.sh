#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 F10V_DATA_ROOT=/workspace/f10v_vision_20261017/data F10V_OUTPUT_ROOT=/workspace/f10v_vision_20261017/runs/prototype_7aed F10V_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
cd /workspace/f10v_vision_20261017/snapshots/ed36aa804/external_connect/vision_tracking_research/collective_motion_v2
while supervisorctl status f10v_train_20261017 | grep -q RUNNING; do sleep 2; done
/workspace/f10v_vision_20261017/.venv/bin/python -m unittest discover -s tests -v
/workspace/f10v_vision_20261017/.venv/bin/python scripts/dependency_audit.py
/workspace/f10v_vision_20261017/.venv/bin/python scripts/evaluate.py freeze
/workspace/f10v_vision_20261017/.venv/bin/python scripts/evaluate.py test
/workspace/f10v_vision_20261017/.venv/bin/python scripts/baselines.py --split final_test --workers 4
/workspace/f10v_vision_20261017/.venv/bin/python scripts/analyze.py
/workspace/f10v_vision_20261017/.venv/bin/python scripts/replay_runtime.py

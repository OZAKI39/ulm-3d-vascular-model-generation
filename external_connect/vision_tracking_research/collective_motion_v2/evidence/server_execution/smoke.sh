#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 F10V_DATA_ROOT=/workspace/f10v_vision_20261017/data F10V_OUTPUT_ROOT=/workspace/f10v_vision_20261017/runs/prototype_7aed
cd /workspace/f10v_vision_20261017/snapshots/7aed42436/external_connect/vision_tracking_research/collective_motion_v2
/workspace/f10v_vision_20261017/.venv/bin/python scripts/audit_environment.py
/workspace/f10v_vision_20261017/.venv/bin/python -m unittest discover -s tests -v
/workspace/f10v_vision_20261017/.venv/bin/python scripts/prepare_data.py --splits train validation calibration final_test --workers 8
/workspace/f10v_vision_20261017/.venv/bin/python scripts/smoke_profile.py

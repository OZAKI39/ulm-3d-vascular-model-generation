#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg
export F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10X_SERVER=/workspace/f10x_piv_20261017
export F10X_OUTPUT=/workspace/f10x_piv_20261017/runs/benchmark_v2
source /workspace/f10x_piv_20261017/.venv/bin/activate
cd /workspace/f10x_piv_20261017/snapshots/3d185b298/external_connect/vision_tracking_research/collective_motion_v4
python -m pip install PySide6-Essentials && python scripts/audit_source.py && python scripts/backend_benchmark.py

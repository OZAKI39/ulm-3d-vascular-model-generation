#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg
export F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
export F10X_SERVER=/workspace/f10x_piv_20261017
export F10X_OUTPUT=/workspace/f10x_piv_20261017/runs/benchmark_checks_v2
source /workspace/f10x_piv_20261017/.venv/bin/activate
cd /workspace/f10x_piv_20261017/snapshots/e05131f38/external_connect/vision_tracking_research/collective_motion_v4
python scripts/backend_benchmark.py --names B3 && python scripts/official_and_consistency.py && /workspace/f10w_vision_20261017/.venv/bin/python scripts/frozen_baseline.py

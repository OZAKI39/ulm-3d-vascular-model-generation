#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg PYTHONPATH=src:.
exec /workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018/.venv/bin/python -m scripts.straight_tube_exp1.f10ab.run analyze_target_confirmation

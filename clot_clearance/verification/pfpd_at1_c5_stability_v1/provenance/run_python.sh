#!/usr/bin/env bash
set -eu
cd /home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONPATH="$PWD/verification/pfpd_at1_petsc_vi_solver_v1/external/python_deps:$PWD${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_CACHE_DIR="$PWD/verification/pfpd_at1_c5_stability_v1/provenance/numba"
export MPLCONFIGDIR="$PWD/verification/pfpd_at1_c5_stability_v1/provenance/mpl"
exec /home/lzy/.cache/codex-dolfinx-env/bin/python "$@"

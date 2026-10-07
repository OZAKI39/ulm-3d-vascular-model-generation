#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export PYTHONPATH=src:.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec .venv/bin/python "$@"

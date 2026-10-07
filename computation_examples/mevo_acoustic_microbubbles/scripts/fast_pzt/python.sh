#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="src:.${PYTHONPATH:+:$PYTHONPATH}"
exec .venv/bin/python "$@"

#!/bin/bash
set -euo pipefail
f10aa_threads=${F10AA_BLAS_THREADS:-1}
case "$f10aa_threads" in 1|2|4) ;; *) echo 'F10AA_BLAS_THREADS must be 1, 2 or 4' >&2; exit 2 ;; esac
export OPENBLAS_NUM_THREADS="$f10aa_threads" OMP_NUM_THREADS="$f10aa_threads" MKL_NUM_THREADS="$f10aa_threads" NUMEXPR_NUM_THREADS=1
export PYTHONPATH="src:.${PYTHONPATH:+:$PYTHONPATH}"
exec .venv/bin/python "$@"

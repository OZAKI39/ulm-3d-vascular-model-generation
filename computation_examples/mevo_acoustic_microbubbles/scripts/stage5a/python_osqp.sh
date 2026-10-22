#!/usr/bin/env bash
set -euo pipefail
stage5a_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
stage5a_python="${MEVO_BASE_PYTHON:-/workspace/mevo_acoustic_microbubbles/.venv/bin/python}"
export PYTHONPATH="${stage5a_root}/.solver_deps/osqp_1p1p3:${stage5a_root}/src:${stage5a_root}${PYTHONPATH:+:${PYTHONPATH}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${stage5a_root}"
exec "${stage5a_python}" "$@"

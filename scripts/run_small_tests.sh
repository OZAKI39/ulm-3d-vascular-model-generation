#!/bin/bash
# Synthetic fixtures and read-only tests only; no BraVa CFD/population batch.
set -euo pipefail
HUMAN_PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HUMAN_PROJECT_ROOT"
TEST_OUTPUT="${1:-validation/small_tests}"
mkdir "$TEST_OUTPUT"
TEST_OUTPUT="$(realpath "$TEST_OUTPUT")"
unset PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CUDA_HOME=/usr/local/cuda-12.6 MAX_JOBS=2
export TORCH_CUDA_ARCH_LIST=8.9
export TORCH_EXTENSIONS_DIR="$HUMAN_PROJECT_ROOT/validation/cache/torch_extensions"
export POPULATION_VALIDATION_DIR="$TEST_OUTPUT/population_artifacts"
export FIXED_VALIDATION_DIR="$TEST_OUTPUT/fixed_window_artifacts"
./.venv/bin/python -I -m pytest tests --junitxml="$TEST_OUTPUT/junit.xml" -v > "$TEST_OUTPUT/pytest.log" 2>&1

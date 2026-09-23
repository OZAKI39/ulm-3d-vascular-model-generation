#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=/workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/repo
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_pipeline events --count 100000 --workers 8 --output /workspace/particle8_2a_20260922/results/admission --provenance /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/provenance.json

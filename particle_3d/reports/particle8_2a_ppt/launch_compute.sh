#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=/workspace/particle8_2a_20260922/83c0cc3ccb5cbc88865ba4e3dcbeb92438f299bd/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/83c0cc3ccb5cbc88865ba4e3dcbeb92438f299bd/repo
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_ppt_data --source /workspace/particle8_2a_20260922/results/trajectories --output /workspace/particle8_2a_20260922/results/ppt --provenance /workspace/particle8_2a_20260922/83c0cc3ccb5cbc88865ba4e3dcbeb92438f299bd/provenance.json --quota 500 --workers 8

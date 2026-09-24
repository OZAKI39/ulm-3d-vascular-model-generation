#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LP_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1
export PYTHONPATH=/workspace/particle8_2a_20260922/7c6343f104ed2c4fec87d17afb2c09c5b8b1a708/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/7c6343f104ed2c4fec87d17afb2c09c5b8b1a708/repo
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals prepare --partial --admission /workspace/particle8_2a_20260922/results/admission --trajectories /workspace/particle8_2a_20260922/results/trajectories --output /workspace/particle8_2a_20260922/results/development_preview
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals render --kind 0 --frames 3 --output /workspace/particle8_2a_20260922/results/development_preview
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals render --kind 1 --frames 6 --output /workspace/particle8_2a_20260922/results/development_preview

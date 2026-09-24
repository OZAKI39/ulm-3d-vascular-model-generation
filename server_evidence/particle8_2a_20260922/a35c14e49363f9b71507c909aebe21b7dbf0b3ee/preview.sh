#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LP_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1
export PYTHONPATH=/workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/repo
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals prepare --partial --admission /workspace/particle8_2a_20260922/results/admission --trajectories /workspace/particle8_2a_20260922/results/trajectories --output /workspace/particle8_2a_20260922/results/development_preview
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals render --kind 0 --frames 3 --output /workspace/particle8_2a_20260922/results/development_preview
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals render --kind 1 --frames 6 --output /workspace/particle8_2a_20260922/results/development_preview
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_analysis admission --development-partial --input /workspace/particle8_2a_20260922/results/admission --report /workspace/particle8_2a_20260922/results/development_preview/review
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_figures --allow-partial --admission /workspace/particle8_2a_20260922/results/admission --report /workspace/particle8_2a_20260922/results/development_preview/review --output /workspace/particle8_2a_20260922/results/development_preview

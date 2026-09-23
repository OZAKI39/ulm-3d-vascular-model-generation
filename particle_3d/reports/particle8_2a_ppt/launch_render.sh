#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=/workspace/particle8_2a_20260922/4a38c85a698c21f4188ed345a7a0b34e947193b6/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/4a38c85a698c21f4188ed345a7a0b34e947193b6/repo
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_ppt_visuals --root /workspace/particle8_2a_20260922/results/ppt
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_ppt_audit --root /workspace/particle8_2a_20260922/results/ppt
/root/particle8_2_runs/env/bin/python -m pytest particle_3d/tests/particle8_2a -q --junitxml=/workspace/particle8_2a_20260922/results/ppt/PPT_REMOTE_TESTS.xml
touch /workspace/particle8_2a_20260922/results/ppt/PPT_COMPLETE

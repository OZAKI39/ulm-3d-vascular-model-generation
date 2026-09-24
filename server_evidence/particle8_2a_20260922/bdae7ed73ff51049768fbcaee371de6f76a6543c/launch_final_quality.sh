#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LP_NUM_THREADS=2 VTK_SMP_MAX_THREADS=2
export PYTHONPATH=/workspace/particle8_2a_20260922/bdae7ed73ff51049768fbcaee371de6f76a6543c/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/bdae7ed73ff51049768fbcaee371de6f76a6543c/repo
while [ ! -f /workspace/particle8_2a_20260922/results/COMPUTE_AND_MEDIA_COMPLETE ]; do sleep 20; done
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals render --output /workspace/particle8_2a_20260922/results/visuals --kind 1
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals audit --output /workspace/particle8_2a_20260922/results/visuals
/root/particle8_2_runs/env/bin/python -m pytest -q particle_3d/tests/particle8_2a --junitxml=/workspace/particle8_2a_20260922/results/review/remote_tests.xml
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_completion_audit --root /workspace/particle8_2a_20260922/results --provenance /workspace/particle8_2a_20260922/bdae7ed73ff51049768fbcaee371de6f76a6543c/provenance.json --admission-provenance /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/provenance.json --trajectory-provenance /workspace/particle8_2a_20260922/fbc4a143212a46dd74069c44dd1207c9b6d29978/provenance.json
touch /workspace/particle8_2a_20260922/results/FINAL_QUALITY_COMPLETE

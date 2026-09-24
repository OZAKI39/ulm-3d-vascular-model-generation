#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LP_NUM_THREADS=2 VTK_SMP_MAX_THREADS=2
export PYTHONPATH=/workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/repo
while [ ! -f /workspace/particle8_2a_20260922/results/admission/ADMISSION_RUN_RECEIPT.json ]; do sleep 20; done
/root/particle8_2_runs/env/bin/python particle_3d/scripts/record_particle82_environment.py --host-provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json --output /workspace/particle8_2a_20260922/results/SERVER_SOFTWARE_ENVIRONMENT.json
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_analysis admission --input /workspace/particle8_2a_20260922/results/admission --report /workspace/particle8_2a_20260922/results/review
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_pipeline benchmark --output /workspace/particle8_2a_20260922/results/scaling_admission_cpu --provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_trajectories --action benchmark --admission /workspace/particle8_2a_20260922/results/admission --output /workspace/particle8_2a_20260922/results/trajectories --provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json
workers=$(/root/particle8_2_runs/env/bin/python -c 'import json;print(json.load(open('"'"'/workspace/particle8_2a_20260922/results/trajectories/TRAJECTORY_SCALING.json'"'"'))["chosen_workers"])')
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_trajectories --action run --admission /workspace/particle8_2a_20260922/results/admission --output /workspace/particle8_2a_20260922/results/trajectories --provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json --workers "$workers"
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_trajectories --action sensitivity --admission /workspace/particle8_2a_20260922/results/admission --output /workspace/particle8_2a_20260922/results/trajectories --provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json --workers "$workers"
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_analysis trajectories --input /workspace/particle8_2a_20260922/results/trajectories --report /workspace/particle8_2a_20260922/results/review
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_sensitivity --admission /workspace/particle8_2a_20260922/results/admission --output /workspace/particle8_2a_20260922/results/sensitivity --provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json --workers "$workers"
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_verify --admission /workspace/particle8_2a_20260922/results/admission --trajectories /workspace/particle8_2a_20260922/results/trajectories --report /workspace/particle8_2a_20260922/results/review --provenance /workspace/particle8_2a_20260922/a35c14e49363f9b71507c909aebe21b7dbf0b3ee/provenance.json
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_figures --admission /workspace/particle8_2a_20260922/results/admission --report /workspace/particle8_2a_20260922/results/review --output /workspace/particle8_2a_20260922/results/visuals
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals prepare --admission /workspace/particle8_2a_20260922/results/admission --trajectories /workspace/particle8_2a_20260922/results/trajectories --output /workspace/particle8_2a_20260922/results/visuals
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals render --output /workspace/particle8_2a_20260922/results/visuals
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_visuals audit --output /workspace/particle8_2a_20260922/results/visuals
/root/particle8_2_runs/env/bin/python -m pytest -q particle_3d/tests/particle8_2a --junitxml=/workspace/particle8_2a_20260922/results/review/remote_tests.xml
touch /workspace/particle8_2a_20260922/results/COMPUTE_AND_MEDIA_COMPLETE

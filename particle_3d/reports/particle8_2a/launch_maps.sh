#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=/workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/repo/particle_3d/src
cd /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/repo
mkdir -p /workspace/particle8_2a_20260922/results/admission
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_pipeline benchmark --output /workspace/particle8_2a_20260922/results/admission --provenance /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/provenance.json
workers=$(/root/particle8_2_runs/env/bin/python -c 'import json;print(json.load(open('"'"'/workspace/particle8_2a_20260922/results/admission/ADMISSION_SCALING.json'"'"'))["chosen_workers"])')
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_pipeline map --map-n 12 --workers "$workers" --output /workspace/particle8_2a_20260922/results/admission --provenance /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/provenance.json
/root/particle8_2_runs/env/bin/python -m particle_3d.particle82a_pipeline map --map-n 24 --workers "$workers" --output /workspace/particle8_2a_20260922/results/admission --provenance /workspace/particle8_2a_20260922/1139128ace574373e3806b6d03866a9b6e49f884/provenance.json
touch /workspace/particle8_2a_20260922/results/admission/MAPS_COMPLETE

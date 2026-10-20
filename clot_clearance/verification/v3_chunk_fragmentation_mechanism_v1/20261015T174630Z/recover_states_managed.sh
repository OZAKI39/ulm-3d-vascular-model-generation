#!/bin/bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd /workspace/y_clot/project/analysis/v3_chunk_fragmentation_mechanism_v1
python -u -B recover_states.py
rc=$?
printf '%s\n' "$rc" > /workspace/y_clot/project/verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z/recover_states_exit_code.txt
exit "$rc"

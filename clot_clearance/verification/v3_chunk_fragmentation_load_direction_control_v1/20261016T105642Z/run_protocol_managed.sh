#!/bin/bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/build:"$LD_LIBRARY_PATH"
cd /workspace/y_clot/project/analysis/v3_chunk_load_direction_control_v1
python -u -B run_protocol.py
rc=$?
printf '%s\n' "$rc" > /workspace/y_clot/project/verification/v3_chunk_fragmentation_load_direction_control_v1/20261016T105642Z/run_protocol_exit_code.txt
exit "$rc"

#!/bin/bash
set -uo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
python -u /workspace/y_clot/project/analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2/run_cases.py --name uniform
uniform_status=$?
echo "UNIFORM_PROCESS_EXIT $uniform_status"
python -u /workspace/y_clot/project/analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2/run_cases.py --name biphasic
biphasic_status=$?
echo "BIPHASIC_PROCESS_EXIT $biphasic_status"
exit $((uniform_status | biphasic_status))

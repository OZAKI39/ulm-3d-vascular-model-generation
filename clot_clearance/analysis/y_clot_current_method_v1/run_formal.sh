#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
exec python -u /workspace/y_clot/project/analysis/y_clot_current_method_v1/run_cases.py --name biphasic

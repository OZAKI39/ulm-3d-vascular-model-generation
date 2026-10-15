#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
TASK_ROOT=/workspace/y_clot/project/verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z
export LD_LIBRARY_PATH="$TASK_ROOT/build:$LD_LIBRARY_PATH"
mpiexec --bind-to core -np 2 python -m analysis.native_state_v1.replay

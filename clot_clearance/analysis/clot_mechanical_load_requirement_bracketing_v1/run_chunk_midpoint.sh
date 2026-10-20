#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/build:${LD_LIBRARY_PATH:-}
while [ ! -f /workspace/y_clot/project/verification/clot_mechanical_load_requirement_bracketing_v1/20261016T202533Z/cases/surface_shear_1000pa_continuation_2/SUMMARY.json ]; do sleep 10; done
exec python -u /workspace/y_clot/project/analysis/clot_mechanical_load_requirement_bracketing_v1/chunk_midpoint.py

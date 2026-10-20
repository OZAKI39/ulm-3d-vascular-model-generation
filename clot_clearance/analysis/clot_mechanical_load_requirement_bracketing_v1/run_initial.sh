#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/build:${LD_LIBRARY_PATH:-}
exec python -u /workspace/y_clot/project/analysis/clot_mechanical_load_requirement_bracketing_v1/campaign.py

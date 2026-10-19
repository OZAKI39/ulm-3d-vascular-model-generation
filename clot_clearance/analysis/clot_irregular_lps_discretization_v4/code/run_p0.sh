#!/usr/bin/env bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
python -u analysis/clot_irregular_lps_discretization_v4/code/p0_prepare.py
export LD_LIBRARY_PATH="/workspace/y_clot/project/verification/clot_irregular_lps_discretization_v4/20261019T221000Z/p0/build:$LD_LIBRARY_PATH"
python -u analysis/clot_irregular_lps_discretization_v4/code/profile_window.py intact
python -u analysis/clot_irregular_lps_discretization_v4/code/profile_window.py damaged

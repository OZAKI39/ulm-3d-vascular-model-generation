#!/bin/bash
set -e
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project/analysis/v3_softening_pilot_v1
exec python -u launch.py pilot

#!/usr/bin/env bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
exec python -u analysis/clot_discretization_current_method_rotation_v1/render.py

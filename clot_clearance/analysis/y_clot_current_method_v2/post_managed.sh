#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
exec python -u /workspace/y_clot/project/analysis/y_clot_current_method_v2/wait_then_post.py

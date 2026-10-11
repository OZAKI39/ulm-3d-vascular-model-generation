#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
exec python -m analysis.lps_pvce_v2.campaign onset

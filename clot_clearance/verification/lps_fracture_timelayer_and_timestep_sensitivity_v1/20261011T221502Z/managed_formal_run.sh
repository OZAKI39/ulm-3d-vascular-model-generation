#!/bin/bash
set -e
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
exec python -u analysis/timelayer_v1/launch.py

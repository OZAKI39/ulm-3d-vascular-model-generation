#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
TASK_DIR=verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z
python analysis/native_state_v1/make_patch.py
python analysis/native_state_v1/compile.py

mpicxx -std=c++17 -O2 -shared -fPIC -DLAMMPS_SMALLBIG -DOMPI_SKIP_MPICXX -I"$TASK_DIR/source/src" -I"$TASK_DIR/source/third_party" analysis/native_state_v1/probe.cpp -L"$TASK_DIR/build" -llammps -o "$TASK_DIR/probe_fixed.so"
export LD_LIBRARY_PATH="/workspace/y_clot/project/$TASK_DIR/build:$LD_LIBRARY_PATH"
python -m analysis.native_state_v1.diagnose sort fixed
mpiexec -np 2 python -m analysis.native_state_v1.diagnose migration fixed
python -m analysis.native_state_v1.diagnose fracture fixed

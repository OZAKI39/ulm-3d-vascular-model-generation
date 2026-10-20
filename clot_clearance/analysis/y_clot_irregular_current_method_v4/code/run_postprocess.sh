#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_irregular_lps_discretization_v4/20261019T221000Z/build:${LD_LIBRARY_PATH:-}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project/analysis/y_clot_irregular_current_method_v4/code
python -u terminal_audit.py
python -u analyze_result.py
python -u render_results.py
python -u finalize.py
python -u validate_delivery.py
python -c "from common import *; dump(O/'POSTPROCESS_STATUS.json',dict(status='COMPLETE',validated=True))"

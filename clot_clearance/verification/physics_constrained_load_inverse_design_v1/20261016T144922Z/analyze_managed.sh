#!/bin/bash
set -euo pipefail
source /workspace/mevo_acoustic_microbubbles/activate.sh
python /workspace/y_clot/project/analysis/physics_constrained_load_inverse_design_v1/mevo_feasibility.py
source /workspace/y_clot/env.sh
exec python /workspace/y_clot/project/analysis/physics_constrained_load_inverse_design_v1/analyze.py

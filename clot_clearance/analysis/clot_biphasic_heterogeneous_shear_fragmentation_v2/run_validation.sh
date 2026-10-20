#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project/analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2
HET_LIB=/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/baseline_capacity/liblammps.so LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/baseline_capacity:${LD_LIBRARY_PATH:-} python -u uniform_probe.py old_uniform
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/build:${LD_LIBRARY_PATH:-}
python -u uniform_probe.py new_uniform
python -u checks.py
python -u path_audit.py
python -u freeze_sources.py

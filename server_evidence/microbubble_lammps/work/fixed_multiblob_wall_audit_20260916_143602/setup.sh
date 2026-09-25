#!/bin/bash
set -euo pipefail
python3 -m venv /workspace/microbubble_lammps/work/fixed_multiblob_wall_audit_20260916_143602/env
/workspace/microbubble_lammps/work/fixed_multiblob_wall_audit_20260916_143602/env/bin/pip install numpy==1.26.4 scipy==1.11.4 numba==0.59.1 h5py==3.11.0 psutil==5.9.8

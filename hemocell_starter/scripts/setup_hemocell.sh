#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec /usr/bin/env -u LD_LIBRARY_PATH -u PYTHONPATH -u PYTHONHOME -u HDF5_DIR \
  PATH=/usr/bin:/bin OMP_NUM_THREADS=1 OMPI_CC=/usr/bin/gcc OMPI_CXX=/usr/bin/g++ \
  /usr/bin/python3 -B "$task_root/scripts/install.py"

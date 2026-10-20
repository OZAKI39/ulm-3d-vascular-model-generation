# Source this file from Bash; only this shell is affected.
_clot_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CLOT_LAMMPS_LIBRARY="${CLOT_LAMMPS_LIBRARY:-${_clot_root}/native/build/liblammps.so}"
export LD_LIBRARY_PATH="$(dirname -- "$CLOT_LAMMPS_LIBRARY")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYVISTA_OFF_SCREEN=true
# Use the active environment's installed dependencies, never a legacy PYTHONPATH.
unset PYTHONPATH
unset _clot_root

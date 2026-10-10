# Source this file to enter the isolated mevo F5 environment.
export MEVO_ROOT=/workspace/mevo_acoustic_microbubbles
export MEVO_PROJECT_DIR="$MEVO_ROOT/repo/computation_examples/mevo_acoustic_microbubbles"
source "$MEVO_ROOT/.venv/bin/activate"
export PYTHONPATH="$MEVO_PROJECT_DIR/src:$MEVO_PROJECT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export WARP_CACHE_PATH="$MEVO_ROOT/.cache/warp"
export XDG_CACHE_HOME="$MEVO_ROOT/.cache"
export MPLCONFIGDIR="$MEVO_ROOT/.cache/matplotlib"
export PYVISTA_OFF_SCREEN=true
cd "$MEVO_PROJECT_DIR"

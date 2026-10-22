#!/bin/bash
set -euo pipefail
app_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4
backend_args=()
if [[ -f "$app_dir/vendor/third_party/openpiv-python-gpu/openpiv_gpu/gpu_process.py" ]]; then
  export F10X_SERVER="$app_dir/vendor"
  backend_args=(--backend B4)
fi
exec "$app_dir/.venv/bin/python" "$app_dir/run_live.py"  "${backend_args[@]}" "$@"

#!/bin/bash
set -euo pipefail
app_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
"$app_dir/.venv/bin/python" -m pip install 'cupy-cuda12x==13.6.0'
mkdir -p "$app_dir/vendor/third_party/openpiv-python-gpu"
tar -xzf "$app_dir/evidence/cupy_runtime_source.tar.gz" -C "$app_dir/vendor/third_party/openpiv-python-gpu"
# A CUDA 12 toolkit/NVRTC and compatible host driver are prerequisites.
# No host driver, camera setting, firewall or system Python is modified.

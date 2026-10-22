#!/bin/bash
set -euo pipefail
F10W_ROOT=/workspace/f10w_vision_20261017
/workspace/f10v_vision_20261017/.venv/bin/python -m venv "$F10W_ROOT/.venv"
cp -a --reflink=auto /workspace/f10v_vision_20261017/.venv/lib/python3.12/site-packages/. "$F10W_ROOT/.venv/lib/python3.12/site-packages/"
uv pip install --python "$F10W_ROOT/.venv/bin/python" "numpy==1.26.4" "POT==0.9.7" pyarrow einops imageio gdown timm hdf5storage
"$F10W_ROOT/.venv/bin/python" -m pip check
"$F10W_ROOT/.venv/bin/python" -m pip freeze > "$F10W_ROOT/setup/requirements.lock"
cd "$F10W_ROOT/third_party"
for repo in LRAFT-ST co-tracker UnLiteFlowNet-PIV; do
 case "$repo" in
 LRAFT-ST) url=https://github.com/wuwuwuas/LRAFT-ST.git;;
 co-tracker) url=https://github.com/facebookresearch/co-tracker.git;;
 UnLiteFlowNet-PIV) url=https://github.com/erizmr/UnLiteFlowNet-PIV.git;;
 esac
 git clone --depth 1 "$url" "$repo"
done

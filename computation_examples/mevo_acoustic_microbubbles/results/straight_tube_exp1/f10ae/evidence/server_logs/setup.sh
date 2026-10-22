#!/bin/bash
set -euo pipefail
cd /workspace/mevo_acoustic_microbubbles/f10ae_target_20261019
uv venv --python /usr/bin/python3.12 .venv
uv pip install --python .venv/bin/python numpy==2.2.6 scipy==1.15.3 h5py==3.16.0 /workspace/mevo_acoustic_microbubbles/migration/wheels/warp_lang-1.18.0+cu12-py3-none-manylinux_2_28_x86_64.whl osqp==1.0.5 matplotlib==3.11.2 pillow==13.0.0 av==19.0.1 pymupdf==1.28.2 pyyaml pytest
uv pip install --python .venv/bin/python torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu
uv pip install --python .venv/bin/python ax-platform==1.3.1 nevergrad==1.0.12
uv pip freeze --python .venv/bin/python > logs/requirements.txt

#!/bin/bash
set -euo pipefail
cd /workspace/mevo_acoustic_microbubbles/f10ae_target_20261019
export UV_HTTP_TIMEOUT=60
uv pip install --python .venv/bin/python torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu
uv pip install --python .venv/bin/python ax-platform==1.3.1 nevergrad==1.0.12
uv pip freeze --python .venv/bin/python > logs/requirements.txt

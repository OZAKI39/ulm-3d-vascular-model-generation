#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export QT_QPA_PLATFORM=xcb
export PYOPENGL_PLATFORM=glx
export PYTHONDONTWRITEBYTECODE=1
exec .venv_annotation/bin/python -B scripts/annotate.py "$@"

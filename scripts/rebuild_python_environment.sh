#!/bin/bash
# New environment only. Does not copy or replace an existing environment.
set -euo pipefail
HUMAN_PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NEW_ENV="${1:-$HUMAN_PROJECT_ROOT/.venv}"
if test -e "$NEW_ENV"; then
    echo 'Environment already exists. Supply a new path; it will not be overwritten.' >&2
    exit 1
fi
unset PYTHONPATH PYTHONHOME
/usr/bin/python3.12 -m venv --copies "$NEW_ENV"
"$NEW_ENV/bin/python" -m pip install --timeout 180 --retries 8 --no-deps \
    'torch==2.6.0+cu126' --index-url https://download.pytorch.org/whl/cu126
# Dependencies use PyPI, avoiding the slower NVIDIA redirect used by the CUDA index.
LOCK_WITHOUT_TORCH=$(mktemp)
trap 'rm -f "$LOCK_WITHOUT_TORCH"' EXIT
"$NEW_ENV/bin/python" - "$HUMAN_PROJECT_ROOT/requirements-lock.txt" "$LOCK_WITHOUT_TORCH" <<'PY'
from pathlib import Path
import sys
lines=[line for line in Path(sys.argv[1]).read_text().splitlines()
       if not line.startswith('torch==') and not line.startswith('-e ')]
Path(sys.argv[2]).write_text('\n'.join(lines)+'\n')
PY
"$NEW_ENV/bin/python" -m pip install --timeout 180 --retries 8 -r "$LOCK_WITHOUT_TORCH"
"$NEW_ENV/bin/python" -m pip install --no-deps --no-build-isolation -e "$HUMAN_PROJECT_ROOT"
"$NEW_ENV/bin/python" -m pip check

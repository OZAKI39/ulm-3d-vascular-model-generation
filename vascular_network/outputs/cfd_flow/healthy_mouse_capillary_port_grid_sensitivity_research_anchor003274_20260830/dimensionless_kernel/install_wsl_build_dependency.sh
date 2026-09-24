#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
printf 'RUN_CWD=%s\n' "$PWD"

python3 -m pip install --user six
python3 -c 'import six; print("SIX_VERSION=" + six.__version__)'


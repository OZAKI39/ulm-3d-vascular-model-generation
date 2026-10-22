#!/usr/bin/env bash
set -euo pipefail
PRINTING_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "$PRINTING_ROOT/config/server_runtime.env" ]]; then
  echo 'Missing config/server_runtime.env: use run_python.sh on the original workstation.' >&2
  exit 1
fi
source "$PRINTING_ROOT/config/server_runtime.env"
exec "$PRINTING_ROOT/run_python.sh" "$@"

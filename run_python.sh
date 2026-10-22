#!/usr/bin/env bash
set -euo pipefail
PRINTING_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$PRINTING_ROOT"
export PYTHONDONTWRITEBYTECODE=1
exec "$PRINTING_ROOT/.venv/bin/python" "$@"

#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$project_dir/microscope_sessions/microscope_20261019_v1/start_annotation.sh" "$@"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"

test -f seeder.lua
mkdir -p mesh

exec /home/lzy/apes-pinned/seeder_official/build/seeder seeder.lua

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"

test -f seeder.lua
test -d ../geometry
test -f ../geometry/geometry_solver_m/wall.stl
test -x /home/lzy/apes-pinned/seeder_official/build/seeder

mkdir -p mesh

printf 'PREFLIGHT_CWD=%s\n' "$PWD"
printf 'SEEDER_LUA_REALPATH=%s\n' "$(realpath seeder.lua)"
printf 'SEEDER_LUA_SHA256=%s\n' "$(sha256sum seeder.lua | awk '{print $1}')"
printf 'WALL_STL_SHA256=%s\n' "$(sha256sum ../geometry/geometry_solver_m/wall.stl | awk '{print $1}')"
printf 'SEEDER_BINARY_SHA256=%s\n' "$(sha256sum /home/lzy/apes-pinned/seeder_official/build/seeder | awk '{print $1}')"
printf 'SEEDER_SOURCE_HEAD=%s\n' "$(git -C /home/lzy/apes-pinned/seeder_official rev-parse HEAD)"

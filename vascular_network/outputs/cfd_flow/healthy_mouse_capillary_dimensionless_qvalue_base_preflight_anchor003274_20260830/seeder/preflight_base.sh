#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
printf 'RUN_CWD=%s\n' "$PWD"

mapfile -t candidates < <(
  find "/home/${USER}/apes-worktrees" -mindepth 1 -maxdepth 1 -type d \
    -name 'seeder_dimensionless_kernel_*' -print | sort
)
[[ ${#candidates[@]} -eq 1 ]]
BINARY="${candidates[0]}/build/seeder"

[[ -f seeder.lua ]]
declare -A expected_geometry=(
  [../geometry/numerical_inlet_plane.stl]='fd071f51e2229971f7b80aa599c78bb8ae3a9d593e06e26460631bcd73b07671'
  [../geometry/geometry_solver_m/outlet_01.stl]='0a204ccea294dc10e1a7298a12e90a4ce14bffa11fc7576ca17d48db88d87593'
  [../geometry/geometry_solver_m/outlet_02.stl]='c884129c40dcc2a1aad6a4d313040be7f7fe3862d376cf2d2dd37488daf457fa'
  [../geometry/geometry_solver_m/outlet_03.stl]='9c51ba9df9d01681180ea7e87edf196ed5ebc179c5a9baaea2c7101ce9306c0a'
  [../geometry/geometry_solver_m/wall.stl]='639ab952e19f60f2f1d1322a4dacbf61d74548eb7b3e37a48834bade7f89bfa5'
)
for geometry in "${!expected_geometry[@]}"; do
  [[ -f "$geometry" ]]
  [[ "$(sha256sum "$geometry" | awk '{print $1}')" == \
    "${expected_geometry[$geometry]}" ]]
done

[[ -x "$BINARY" ]]
[[ "$(sha256sum "$BINARY" | awk '{print $1}')" == \
  'd7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed' ]]
[[ "$(git -C "${candidates[0]}" rev-parse HEAD)" == \
  '667109df6fafdcb39f4409e3f5d90f04d75cd33c' ]]
[[ "$(git -C "${candidates[0]}/tem" rev-parse HEAD)" == \
  '53f273dbb8e9dcbe7feeb3d9831a35f5ae3cd72c' ]]

mkdir -p mesh
probe_file="mesh/.write_probe_$$"
: > "$probe_file"
rm -- "$probe_file"

printf 'CONFIG=%s\n' "$SCRIPT_DIR/seeder.lua"
printf 'CONFIG_SHA256=%s\n' "$(sha256sum seeder.lua | awk '{print $1}')"
printf 'BINARY=%s\n' "$BINARY"
printf 'BINARY_SHA256=%s\n' "$(sha256sum "$BINARY" | awk '{print $1}')"
printf 'SOURCE_HEAD=%s\n' "$(git -C "${candidates[0]}" rev-parse HEAD)"
printf 'OUTPUT_DIR=%s\n' "$SCRIPT_DIR/mesh"
printf 'PREFLIGHT=PASS\n'

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
[[ -f geometry/wall.stl ]]
[[ -x "$BINARY" ]]
[[ "$(sha256sum geometry/wall.stl | awk '{print $1}')" == \
  '8e99cbd7564ee36d7ce751538e81f2a659a971df68038ebbfdb2603ab67ed069' ]]
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
printf 'GEOMETRY=%s\n' "$SCRIPT_DIR/geometry/wall.stl"
printf 'BINARY=%s\n' "$BINARY"
printf 'BINARY_SHA256=%s\n' "$(sha256sum "$BINARY" | awk '{print $1}')"
printf 'SOURCE_HEAD=%s\n' "$(git -C "${candidates[0]}" rev-parse HEAD)"
printf 'OUTPUT_DIR=%s\n' "$SCRIPT_DIR/mesh"
printf 'PREFLIGHT=PASS\n'

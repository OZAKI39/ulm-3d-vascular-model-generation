#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
printf 'RUN_CWD=%s\n' "$PWD"

TARGET_ROOT="${HOME}/apes-worktrees/seeder_dimensionless_kernel_20260830"
test -d "$TARGET_ROOT"

git -C "$TARGET_ROOT" submodule update --init --recursive
test -f "$TARGET_ROOT/tem/source/shapes/tem_line_module.fpp"
test -f "$TARGET_ROOT/sdr/source/seeder.f90"

printf 'SEEDER_HEAD=%s\n' "$(git -C "$TARGET_ROOT" rev-parse HEAD)"
printf 'TREELM_HEAD=%s\n' "$(git -C "$TARGET_ROOT/tem" rev-parse HEAD)"
printf 'SDR_HEAD=%s\n' "$(git -C "$TARGET_ROOT/sdr" rev-parse HEAD)"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
printf 'RUN_CWD=%s\n' "$PWD"

PINNED_ROOT="${HOME}/apes-pinned/seeder_official"
TARGET_ROOT="${HOME}/apes-worktrees/seeder_dimensionless_kernel_20260830"
test -d "$PINNED_ROOT"
test -d "$TARGET_ROOT"

for module in aotus polynomials bin tem sdr; do
  test -d "$PINNED_ROOT/$module"
  git -C "$TARGET_ROOT" config "submodule.${module}.url" "$PINNED_ROOT/$module"
done

git -c protocol.file.allow=always -C "$TARGET_ROOT" submodule update --init
test -f "$TARGET_ROOT/tem/source/shapes/tem_line_module.fpp"
test -f "$TARGET_ROOT/sdr/source/seeder.f90"

printf 'SEEDER_HEAD=%s\n' "$(git -C "$TARGET_ROOT" rev-parse HEAD)"
printf 'TREELM_HEAD=%s\n' "$(git -C "$TARGET_ROOT/tem" rev-parse HEAD)"
printf 'SDR_HEAD=%s\n' "$(git -C "$TARGET_ROOT/sdr" rev-parse HEAD)"

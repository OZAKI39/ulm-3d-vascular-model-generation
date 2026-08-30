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
BASE_SHA='667109df6fafdcb39f4409e3f5d90f04d75cd33c'

test -d "$PINNED_ROOT"
test "$(git -C "$PINNED_ROOT" rev-parse HEAD)" = "$BASE_SHA"
test ! -e "$TARGET_ROOT"

git -C "$PINNED_ROOT" worktree add --detach "$TARGET_ROOT" "$BASE_SHA"
git -C "$TARGET_ROOT" submodule update --init --recursive

printf 'TARGET_ROOT=%s\n' "$TARGET_ROOT"
printf 'SEEDER_HEAD=%s\n' "$(git -C "$TARGET_ROOT" rev-parse HEAD)"
printf 'TREELM_HEAD=%s\n' "$(git -C "$TARGET_ROOT/tem" rev-parse HEAD)"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
printf 'RUN_CWD=%s\n' "$PWD"

PROJECT_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
PATCH_FILE="$PROJECT_ROOT/patches/seeder/dimensionless_ray_triangle.patch"
WORKTREE_ROOT="${SEEDER_DIMENSIONLESS_WORKTREE:-}"

if [[ -z "$WORKTREE_ROOT" ]]; then
  mapfile -t worktree_candidates < <(
    find "/home/${USER}/apes-worktrees" -mindepth 1 -maxdepth 1 -type d \
      -name 'seeder_dimensionless_kernel_*' -print | sort
  )
  if [[ ${#worktree_candidates[@]} -ne 1 ]]; then
    printf 'Expected one dimensionless Seeder worktree; found %d\n' \
      "${#worktree_candidates[@]}" >&2
    exit 2
  fi
  WORKTREE_ROOT="${worktree_candidates[0]}"
fi

SOURCE_FILE="$WORKTREE_ROOT/tem/source/shapes/tem_line_module.fpp"
[[ -f "$PATCH_FILE" ]]
[[ -f "$SOURCE_FILE" ]]

EXPECTED_SEEDER_SHA='667109df6fafdcb39f4409e3f5d90f04d75cd33c'
EXPECTED_TREELM_SHA='53f273dbb8e9dcbe7feeb3d9831a35f5ae3cd72c'
ACTUAL_SEEDER_SHA="$(git -C "$WORKTREE_ROOT" rev-parse HEAD)"
ACTUAL_TREELM_SHA="$(git -C "$WORKTREE_ROOT/tem" rev-parse HEAD)"
[[ "$ACTUAL_SEEDER_SHA" == "$EXPECTED_SEEDER_SHA" ]]
[[ "$ACTUAL_TREELM_SHA" == "$EXPECTED_TREELM_SHA" ]]

# A reverse check proves that the complete recorded patch, and no hand-edited
# approximation of it, is present in the external TreElm source tree.
git -C "$WORKTREE_ROOT/tem" apply --check --reverse "$PATCH_FILE"

if [[ -x "/home/${USER}/.local/bin/mpif90" ]]; then
  MPI_FC="/home/${USER}/.local/bin/mpif90"
else
  MPI_FC="$(command -v mpif90)"
fi

mkdir -p "$SCRIPT_DIR/build"
BUILD_LOG="$SCRIPT_DIR/build/seeder_dimensionless_build.log"
{
  printf 'PROJECT_ROOT=%s\n' "$PROJECT_ROOT"
  printf 'WORKTREE_ROOT=%s\n' "$WORKTREE_ROOT"
  printf 'SEEDER_BASE_SHA=%s\n' "$ACTUAL_SEEDER_SHA"
  printf 'TREELM_BASE_SHA=%s\n' "$ACTUAL_TREELM_SHA"
  printf 'PATCH_SHA256=%s\n' "$(sha256sum "$PATCH_FILE" | awk '{print $1}')"
  printf 'FC=%s\n' "$MPI_FC"
  cd "$WORKTREE_ROOT"
  FC="$MPI_FC" bin/waf configure build
  [[ -x "$WORKTREE_ROOT/build/seeder" ]]
  printf 'BINARY_PATH=%s\n' "$WORKTREE_ROOT/build/seeder"
  printf 'BINARY_SHA256=%s\n' \
    "$(sha256sum "$WORKTREE_ROOT/build/seeder" | awk '{print $1}')"
} 2>&1 | tee "$BUILD_LOG"

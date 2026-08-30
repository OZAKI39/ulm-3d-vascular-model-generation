#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
SOURCE_ROOT='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300'
MUS_ROOT="$SOURCE_ROOT/mus"
TEM_ROOT="$SOURCE_ROOT/tem"
BINARY="$SOURCE_ROOT/build/musubi_adaptive_flux"
OUT_DIR="$SCRIPT_DIR/source_trace"

[[ -d "$MUS_ROOT/.git" || -f "$MUS_ROOT/.git" ]]
[[ -d "$TEM_ROOT/.git" || -f "$TEM_ROOT/.git" ]]
[[ -x "$BINARY" ]]
mkdir -p "$OUT_DIR"

{
  printf 'SOURCE_ROOT=%s\n' "$SOURCE_ROOT"
  printf 'OUTER_HEAD=%s\n' "$(git -C "$SOURCE_ROOT" rev-parse HEAD)"
  printf 'MUS_HEAD=%s\n' "$(git -C "$MUS_ROOT" rev-parse HEAD)"
  printf 'TEM_HEAD=%s\n' "$(git -C "$TEM_ROOT" rev-parse HEAD)"
  printf 'BINARY=%s\n' "$BINARY"
  printf 'BINARY_SHA256=%s\n' "$(sha256sum "$BINARY" | awk '{print $1}')"
} > "$OUT_DIR/revisions.log"

git -C "$MUS_ROOT" grep -n -E \
  'glob_source|force_order|body_force|applySrc_force|derive_force|forceField|velocity.*force|physics%fac.*body_force|mus_calcOmegaFromVisc|cs2inv.*visc' \
  -- source > "$OUT_DIR/force_git_grep.log"

{
  git -C "$MUS_ROOT" log --oneline -S'body_force' -- source/mus_physics_module.f90
  git -C "$MUS_ROOT" log --oneline -S'force_order' -- source
  git -C "$MUS_ROOT" log --oneline -G'applySrc_force|derive_force' -- source
} > "$OUT_DIR/force_git_history.log"

git -C "$MUS_ROOT" grep -n -E \
  'wall_libb|qVal|iMeshDir|stencil%map|invDir|cxDirInv|PULL|FETCH|bouzidi|fIn|fOut|fNgh' \
  -- source > "$OUT_DIR/wall_git_grep.log"

{
  git -C "$MUS_ROOT" log --oneline -S'wall_libb' -- source
  git -C "$MUS_ROOT" log --oneline -S'iMeshDir = stencil%map' -- source
  git -C "$MUS_ROOT" log --oneline -G'invDir|cxDirInv|bouzidi' -- source
  git -C "$TEM_ROOT" log --oneline -G'qVal|property' -- source
} > "$OUT_DIR/wall_git_history.log"

{
  sed -n '520,545p' "$MUS_ROOT/source/mus_physics_module.f90"
  sed -n '3108,3245p' "$MUS_ROOT/source/derived/mus_derQuan_module.fpp"
  sed -n '1150,1220p' "$MUS_ROOT/source/derived/mus_auxFieldVar_module.fpp"
  sed -n '350,368p' "$MUS_ROOT/source/mus_relaxationParam_module.f90"
} > "$OUT_DIR/force_source_excerpts.log"

{
  sed -n '320,405p' "$MUS_ROOT/source/bc/mus_bc_fluid_wall_module.fpp"
  sed -n '1750,1839p' "$MUS_ROOT/source/bc/mus_bc_header_module.fpp"
  sed -n '1345,1383p' "$MUS_ROOT/source/bc/mus_bc_general_module.fpp"
  sed -n '2300,2350p' "$MUS_ROOT/source/mus_construction_module.fpp"
  sed -n '1,180p' "$MUS_ROOT/source/header/lbm_macros.inc"
} > "$OUT_DIR/wall_source_excerpts.log"

for required in revisions.log force_git_grep.log force_git_history.log \
  wall_git_grep.log wall_git_history.log force_source_excerpts.log \
  wall_source_excerpts.log; do
  [[ -s "$OUT_DIR/$required" ]]
done

printf 'SOURCE_TRACE_SEMANTIC_SUCCESS=PASS\n'

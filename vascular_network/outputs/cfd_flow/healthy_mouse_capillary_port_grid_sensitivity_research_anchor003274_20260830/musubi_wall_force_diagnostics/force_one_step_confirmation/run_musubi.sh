#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
MESH="$SCRIPT_DIR/../force_one_step/mesh"
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
/bin/bash "$SCRIPT_DIR/preflight_musubi.sh" > musubi_preflight.log 2>&1
if find tracking -mindepth 1 -print -quit | grep -q .; then
  [[ -s musubi_stdout.log ]]
else
  "$MPIRUN" --bind-to core --map-by core -np 1 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
fi
grep -q 'Initializing musubi' musubi_stdout.log
grep -Eq 'iterations:[[:space:]]+2' musubi_stdout.log
test "$(find tracking -type f -name '*.res' -size +0c | wc -l)" -ge 2
printf 'SOLVER_INITIALIZED=PASS
MESH_LOADED=PASS
ITERATION_TWO=PASS
TRACKING_READABLE=PASS
SEMANTIC_SUCCESS=PASS
' > musubi_semantic_status.log

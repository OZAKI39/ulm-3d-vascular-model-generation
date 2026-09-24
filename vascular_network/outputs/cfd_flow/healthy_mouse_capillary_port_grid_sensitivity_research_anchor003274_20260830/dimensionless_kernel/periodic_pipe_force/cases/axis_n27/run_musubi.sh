#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$SCRIPT_DIR"
SEEDER='/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
FREEZE="$SCRIPT_DIR/../../../../qc/seeder_geometry_freeze.json"
/bin/bash "$SCRIPT_DIR/preflight_musubi.sh" > musubi_preflight.log 2>&1
if find tracking -mindepth 1 -print -quit | grep -q .; then
  printf 'Refusing to overwrite non-empty tracking directory
' >&2; exit 3
fi
"$MPIRUN" --bind-to core --map-by core -np 2 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
grep -q 'Initializing simulation' musubi_stdout.log
test "$(find tracking -type f -name '*.res' -size +0c | wc -l)" -ge 3
printf 'SEMANTIC_SUCCESS=PASS
' > musubi_semantic_status.log

#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
SEEDER='/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
/bin/bash "$SCRIPT_DIR/preflight_seeder.sh" > seeder_preflight.log 2>&1
if find mesh -mindepth 1 -print -quit | grep -q .; then
  [[ -s mesh/header.lua && -s mesh/elemlist.lsb && -s seeder_stdout.log ]]
else
  "$SEEDER" seeder.lua > seeder_stdout.log 2> seeder_stderr.log
fi
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "mesh/$f" ]]; done
grep -q 'Done with Seeder' seeder_stdout.log
grep -q 'nElems = 32' mesh/header.lua
printf 'MESH_LOADED=PASS\nCELL_COUNT=32\nQVAL_ARTIFACT=PASS\nSEMANTIC_SUCCESS=PASS\n' > seeder_semantic_status.log

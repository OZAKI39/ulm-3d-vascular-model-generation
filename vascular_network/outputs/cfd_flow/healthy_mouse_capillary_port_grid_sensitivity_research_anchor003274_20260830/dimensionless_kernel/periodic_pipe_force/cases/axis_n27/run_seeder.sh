#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$SCRIPT_DIR"
SEEDER='/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
FREEZE="$SCRIPT_DIR/../../../../qc/seeder_geometry_freeze.json"
/bin/bash "$SCRIPT_DIR/preflight_seeder.sh" > seeder_preflight.log 2>&1
if find mesh -mindepth 1 -print -quit | grep -q .; then
  printf 'Refusing to overwrite non-empty mesh directory
' >&2; exit 3
fi
"$SEEDER" seeder.lua > seeder_stdout.log 2> seeder_stderr.log
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "mesh/$f" ]]; done
grep -q 'Seeder created mesh successfully!' seeder_stdout.log
printf 'SEMANTIC_SUCCESS=PASS
' > seeder_semantic_status.log

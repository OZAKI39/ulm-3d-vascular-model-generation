#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$SCRIPT_DIR"
SEEDER='/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
FREEZE="$SCRIPT_DIR/../../../../qc/seeder_geometry_freeze.json"
[[ -f musubi.lua && -d mesh && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ -f "$FREEZE" ]] && grep -q 'CFD_FLOW_SEEDER_GEOMETRY_KERNEL_VALIDATED' "$FREEZE"
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "mesh/$f" ]]; done
mkdir -p tracking
probe="tracking/.write_probe_$$"; : > "$probe"; rm -- "$probe"
printf 'CONFIG=%s\nBINARY_SHA256=%s\nPREFLIGHT=PASS\n' "$SCRIPT_DIR/musubi.lua" 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588'

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
[[ -f musubi.lua && -d "$MESH" && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
for f in header.lua elemlist.lsb; do [[ -s "$MESH/$f" ]]; done
grep -q 'nElems = 64' "$MESH/header.lua"
mkdir -p tracking
probe="tracking/.write_probe_$$"; : > "$probe"; rm -- "$probe"
printf 'CONFIG=%s\nMESH=%s\nBINARY_SHA256=%s\nPREFLIGHT=PASS\n' "$SCRIPT_DIR/musubi.lua" "$MESH" 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588'

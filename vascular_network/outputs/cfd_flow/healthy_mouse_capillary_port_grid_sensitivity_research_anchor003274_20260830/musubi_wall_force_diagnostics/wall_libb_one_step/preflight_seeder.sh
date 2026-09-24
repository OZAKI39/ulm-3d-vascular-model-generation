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
[[ -f seeder.lua && -s geometry/wall.stl && -x "$SEEDER" ]]
[[ "$(sha256sum "$SEEDER" | awk '{print $1}')" == 'd7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed' ]]
mkdir -p mesh
probe="mesh/.write_probe_$$"; : > "$probe"; rm -- "$probe"
printf 'CONFIG=%s\nBINARY_SHA256=%s\nPREFLIGHT=PASS\n' "$SCRIPT_DIR/seeder.lua" 'd7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed'

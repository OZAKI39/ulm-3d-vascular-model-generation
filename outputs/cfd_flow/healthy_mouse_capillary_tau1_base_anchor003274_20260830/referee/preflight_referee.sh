#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
MESH='/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
[[ -f musubi.lua && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
grep -q 'nElems = 182320' "$MESH/header.lua"
[[ "$(sha256sum "$MESH/elemlist.lsb" | awk '{print $1}')" == 'f7d7b1d55273c78c336ac04e39bc018dd9ebb470a9f29ce833ff01711de8c386' ]]
[[ "$(sha256sum "$MESH/bnd.lsb" | awk '{print $1}')" == '520d7dd1e4a46a45f9b1218a5807cfd89d6f054e0a247872362b130ff6bcfe69' ]]
[[ "$(sha256sum "$MESH/qval.lsb" | awk '{print $1}')" == '35884406b5f0111cd4ab471f7b08ac3df00e478d3458a57636d1bd8921cb0fe6' ]]
mkdir -p tracking restart
printf 'CELL_COUNT=182320
BINARY_SHA256=e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588
PREFLIGHT=PASS
'

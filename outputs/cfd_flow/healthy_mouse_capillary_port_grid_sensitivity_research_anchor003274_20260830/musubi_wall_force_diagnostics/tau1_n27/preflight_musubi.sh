#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
MESH='/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_port_grid_sensitivity_research_anchor003274_20260830/dimensionless_kernel/periodic_pipe_force/cases/axis_n27/mesh'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
[[ -f musubi.lua && -d "$MESH" && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
grep -q 'nElems = 90720' "$MESH/header.lua"
mkdir -p tracking
probe="tracking/.write_probe_$$"; : > "$probe"; rm -- "$probe"
printf 'CONFIG=%s
MESH=%s
BINARY_SHA256=%s
CELL_COUNT=90720
PREFLIGHT=PASS
' "$SCRIPT_DIR/musubi.lua" "$MESH" 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588'

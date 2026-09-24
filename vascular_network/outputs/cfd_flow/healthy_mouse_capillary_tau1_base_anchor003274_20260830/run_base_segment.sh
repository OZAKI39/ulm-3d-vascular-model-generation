#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
SEGMENT_DIR="${1:?segment directory required}"
MESH='/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
[[ -d "$SEGMENT_DIR" && -s "$SEGMENT_DIR/musubi.lua" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
grep -q 'nElems = 182320' "$MESH/header.lua"
mkdir -p '/home/lzy/u3da/tau1_base_20260830/restart' "$SEGMENT_DIR/tracking/p" "$SEGMENT_DIR/tracking/u"
cd "$SEGMENT_DIR"
set +e
"$MPIRUN" --bind-to core --map-by core --report-bindings -np 4 "$MUSUBI" musubi.lua 2> musubi_stderr.log | awk -v interval=59875 '
/ADAPTIVE_FLUX_PRESSURE/ {
  last=$0
  if (match($0,/iter=[0-9]+/)) { value=substr($0,RSTART+5,RLENGTH-5)+0; if (value % interval == 0) { print $0; fflush() } }
  next
}
{ print $0 }
END { if (last != "") print last }
' > musubi_stdout.log
rc=${PIPESTATUS[0]}
set -e
[[ "$rc" -eq 0 ]]
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Got a mesh with following properties' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE' musubi_stdout.log
grep -q 'SUCCESSFUL run' musubi_stdout.log
[[ -s '/home/lzy/u3da/tau1_base_20260830/restart/tau1_base_lastHeader.lua' ]]
printf 'SEGMENT_SEMANTIC_SUCCESS=PASS
' > semantic_status.log

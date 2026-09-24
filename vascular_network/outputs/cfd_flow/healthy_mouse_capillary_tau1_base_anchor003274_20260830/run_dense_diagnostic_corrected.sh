#!/usr/bin/env bash
set -euo pipefail
DENSE_DIR="${1:?dense diagnostic directory required}"
MESH='/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh'
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
[[ -d "$DENSE_DIR" && -s "$DENSE_DIR/musubi.lua" ]]
[[ ! -e "$DENSE_DIR/diagnostic_launched" ]]
[[ "$(sha256sum '/home/lzy/u3da/tau1_base_20260830/restart/tau1_base_6.357E-03.lsb' | awk '{print $1}')" == '3d54f3970b4120896c214155811d7cd1b594e3efd172f80b5dc5e7d0fef279e2' ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
mkdir -p "$DENSE_DIR/restart" "$DENSE_DIR/tracking/p" "$DENSE_DIR/tracking/u"
touch "$DENSE_DIR/diagnostic_launched"
cd "$DENSE_DIR"
"$MPIRUN" --bind-to core --map-by core --report-bindings -np 4 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Got a mesh with following properties' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'Found BC wall of kind wall_libb' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE' musubi_stdout.log
grep -q 'SUCCESSFUL run' musubi_stdout.log
printf 'DENSE_DIAGNOSTIC_SEMANTIC_SUCCESS=PASS
' > semantic_status.log

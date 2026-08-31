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
/bin/bash "$SCRIPT_DIR/preflight_referee.sh" > preflight.log 2>&1
if find restart -maxdepth 1 -type f -name '*.lsb' -size +0c | grep -q .; then
  [[ -s musubi_stdout.log ]]
else
  "$MPIRUN" --bind-to core --map-by core --report-bindings -np 4 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
fi
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Got a mesh with following properties' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'Found BC wall of kind wall_libb' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE iter=1' musubi_stdout.log
grep -Eq 'iterations:[[:space:]]+1' musubi_stdout.log
test "$(find restart -maxdepth 1 -type f -iname '*header*.lua' -size +0c | wc -l)" -ge 2
test "$(find restart -maxdepth 1 -type f -name '*.lsb' -size +0c | wc -l)" -ge 2
printf 'MESH_LOADED=PASS
CONTINUOUS_Q_LOADED=PASS
ITERATION_ONE=PASS
RESTART_ZERO_ONE=PASS
SEMANTIC_SUCCESS=PASS
' > semantic_status.log

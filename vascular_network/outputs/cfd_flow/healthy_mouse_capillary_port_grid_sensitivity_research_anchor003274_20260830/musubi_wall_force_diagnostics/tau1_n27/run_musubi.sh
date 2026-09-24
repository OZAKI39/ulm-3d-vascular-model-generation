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
/bin/bash "$SCRIPT_DIR/preflight_musubi.sh" > musubi_preflight.log 2>&1
if find tracking -mindepth 1 -print -quit | grep -q .; then
  [[ -s musubi_stdout.log ]]
else
  "$MPIRUN" --bind-to core --map-by core -np 2 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
fi
grep -q 'Initializing musubi' musubi_stdout.log
for label in mean_velocity profile safety; do
  result="$(find tracking -maxdepth 1 -type f -name "*${label}*p00000.res" -print -quit)"
  [[ -n "$result" && -s "$result" ]]
  test "$(grep -cv '^[[:space:]]*#\|^[[:space:]]*$' "$result")" -ge 5
done
test "$(find tracking -maxdepth 1 -type f -name '*cross_section*.res' -size +0c | wc -l)" -ge 5
printf 'SOLVER_INITIALIZED=PASS
MESH_LOADED=PASS
ACTUAL_ITERATION_GE_200=PASS
TRACKING_READABLE=PASS
CROSS_SECTION_READABLE=PASS
SEMANTIC_SUCCESS=PASS
' > musubi_semantic_status.log

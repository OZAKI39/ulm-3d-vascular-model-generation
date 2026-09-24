#!/usr/bin/env bash
set -euo pipefail
SEGMENT_DIR="${1:?segment directory required}"
cd "$SEGMENT_DIR"
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
[[ -s musubi.lua && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
mkdir -p tracking '/home/lzy/u3da/tau1_reference_scaled_base_20260901/restart'
set +e
"$MPIRUN" --bind-to core --map-by core --report-bindings -np 4   "$MUSUBI" musubi.lua 2> musubi_stderr.log | awk -v interval=119751 '
/ADAPTIVE_FLUX_PRESSURE/ {
  last=$0
  if (match($0,/iter=[0-9]+/)) {
    value=substr($0,RSTART+5,RLENGTH-5)+0
    if (value % interval == 0) { print $0; fflush() }
  }
  next
}
{ print $0; fflush() }
END { if (last != "") print last }
' > musubi_stdout.log
rc=${PIPESTATUS[0]}
set -e
[[ "$rc" -eq 0 ]]
grep -q 'SUCCESSFUL run' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE' musubi_stdout.log
printf 'SEGMENT_SEMANTIC_SUCCESS=PASS
' > semantic_status.log

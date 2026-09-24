#!/usr/bin/env bash
set -euo pipefail
SEGMENT_DIR="${1:?segment directory required}"
cd "$SEGMENT_DIR"
mkdir -p tracking
set +e
'/home/lzy/.local/bin/mpirun' --bind-to core --map-by core --report-bindings -np 4 \
  '/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux' musubi.lua 2> musubi_stderr.log | awk -v interval=70859 '
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

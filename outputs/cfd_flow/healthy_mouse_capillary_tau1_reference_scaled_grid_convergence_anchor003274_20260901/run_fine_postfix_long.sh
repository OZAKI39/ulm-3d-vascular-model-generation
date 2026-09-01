#!/usr/bin/env bash
set -euo pipefail

segment=/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/segments/segment_0005000_to_1011895
status=/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/long_run_status.txt
printf '%s\n' "$$" > /home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/long_run_pid.txt
printf 'RUNNING\n' > "$status"
cd "$segment"
set +e
/usr/bin/time -v /home/lzy/.local/bin/mpirun --bind-to core --map-by core --report-bindings -np 4 \
  /home/lzy/apes-worktrees/musubi_adaptive_target_fixed_20260901/build/musubi musubi.lua \
  2> musubi_stderr.log | awk -v interval=202379 '
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
if [[ "$rc" -eq 0 ]] && grep -q 'SUCCESSFUL run' musubi_stdout.log; then
  printf 'PASS\n' > "$status"
  exit 0
fi
printf 'FAIL rc=%s\n' "$rc" > "$status"
exit "$rc"

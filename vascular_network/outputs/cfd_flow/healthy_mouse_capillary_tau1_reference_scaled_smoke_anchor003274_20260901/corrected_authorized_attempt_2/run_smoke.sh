#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$SCRIPT_DIR"
MUSUBI='/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/build/musubi_adaptive_flux'
MPIRUN='/home/lzy/.local/bin/mpirun'
[[ -x "$MUSUBI" && -x "$MPIRUN" && -s musubi.lua ]]
[[ "$(sha256sum "$MUSUBI" | awk '{print $1}')" == 'e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588' ]]
! grep -Eq 'restart[[:space:]]*=.*read|read[[:space:]]*=' musubi.lua
if find restart -maxdepth 1 -type f -size +0c 2>/dev/null | grep -q .; then
  echo 'fresh smoke restart directory is not empty' >&2
  exit 9
fi
mkdir -p restart tracking/state_000000 tracking/state_000001 tracking/state_000010 tracking/state_000100 tracking/state_001000 tracking/state_005000
"$MPIRUN" --bind-to core --map-by core --report-bindings -np 4 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'Found BC wall of kind wall_libb' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE iter=5000' musubi_stdout.log
grep -q 'SUCCESSFUL run' musubi_stdout.log
printf 'FRESH_INITIALIZATION=PASS
ITERATIONS=5000
SEMANTIC_SUCCESS=PASS
' > semantic_status.log

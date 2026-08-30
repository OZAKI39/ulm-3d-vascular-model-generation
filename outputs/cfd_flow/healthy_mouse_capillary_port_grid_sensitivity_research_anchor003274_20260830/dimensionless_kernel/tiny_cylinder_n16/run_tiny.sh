#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
)"

cd "$SCRIPT_DIR"
printf 'RUN_CWD=%s\n' "$PWD"

/bin/bash "$SCRIPT_DIR/preflight_tiny.sh" > launcher_preflight.log 2>&1
mapfile -t candidates < <(
  find "/home/${USER}/apes-worktrees" -mindepth 1 -maxdepth 1 -type d \
    -name 'seeder_dimensionless_kernel_*' -print | sort
)
BINARY="${candidates[0]}/build/seeder"

if find mesh -mindepth 1 -print -quit | grep -q .; then
  printf 'Refusing to overwrite a non-empty mesh directory\n' >&2
  exit 3
fi

started="$(date +%s)"
"$BINARY" seeder.lua > seeder_stdout.log 2> seeder_stderr.log
finished="$(date +%s)"

for required in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do
  [[ -s "mesh/$required" ]]
done
grep -q 'Seeder created mesh successfully!' seeder_stdout.log
printf 'SEEDER_RUNTIME_SECONDS=%s\n' "$((finished - started))" > runtime.log
printf 'SEMANTIC_SUCCESS=PASS\n' >> runtime.log

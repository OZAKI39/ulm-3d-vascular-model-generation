#!/usr/bin/env bash
# Source from the desired immutable snapshot or local F10AF worktree.
F10AF_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "$F10AF_CODE_ROOT/.venv/bin/activate"
export PYTHONPATH="$F10AF_CODE_ROOT/src:$F10AF_CODE_ROOT"
export MPLBACKEND=Agg

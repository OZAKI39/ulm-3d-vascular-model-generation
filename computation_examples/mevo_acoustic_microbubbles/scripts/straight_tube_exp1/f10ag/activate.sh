#!/bin/bash
F10AG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "$F10AG_ROOT/.venv/bin/activate"
export PYTHONPATH="$F10AG_ROOT/src:$F10AG_ROOT" MPLBACKEND=Agg

#!/usr/bin/env bash
set -euo pipefail
exec ssh -i /home/lzy/.ssh/vast_step3b_ed25519 -p 15283 -o IdentitiesOnly=yes -o IPQoS=none root@95.3.33.46 "$@"

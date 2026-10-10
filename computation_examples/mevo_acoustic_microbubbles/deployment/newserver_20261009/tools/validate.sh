#!/bin/bash
set -euo pipefail
source /workspace/mevo_acoustic_microbubbles/activate.sh
timeout --signal=TERM --kill-after=30s 1200s python -u "$MEVO_ROOT/migration/tools/verify_migration.py" 2>&1 | tee "$MEVO_ROOT/migration/logs/validation.log"

#!/bin/bash
set -euo pipefail
source /workspace/mevo_acoustic_microbubbles/activate-render.sh
timeout --signal=TERM --kill-after=30s 300s python -u "$MEVO_ROOT/migration/tools/verify_migration.py" --phase render_video 2>&1 | tee "$MEVO_ROOT/migration/logs/render-validation.log"

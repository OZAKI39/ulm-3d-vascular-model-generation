#!/bin/bash
set -euo pipefail
/venv/main/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/gpu_validate.py mesh
/root/particle8_2_runs/env/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/test_native.py
/root/particle8_2_runs/env/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/verify_native_bytes.py
/root/particle8_2_runs/env/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/production_entry.py
/venv/main/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/gpu_validate.py tracks
/root/particle8_2_runs/env/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/summarize.py

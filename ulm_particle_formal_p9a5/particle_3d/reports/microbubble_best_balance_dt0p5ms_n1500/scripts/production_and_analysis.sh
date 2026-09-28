#!/bin/bash
set -euo pipefail
/root/particle8_2_runs/env/bin/python /workspace/microbubble_best_balance_dt0p5ms_n1500_20260928/scripts/production_entry.py
/venv/main/bin/python /workspace/microbubble_best_balance_dt0p5ms_n1500_20260928/scripts/gpu_validate.py tracks
/root/particle8_2_runs/env/bin/python /workspace/microbubble_best_balance_dt0p5ms_n1500_20260928/scripts/summarize.py

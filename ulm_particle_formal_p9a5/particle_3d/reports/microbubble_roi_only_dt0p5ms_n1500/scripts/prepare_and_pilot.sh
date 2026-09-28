#!/bin/bash
set -euo pipefail
/root/particle8_2_runs/env/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/campaign.py prepare
/root/particle8_2_runs/env/bin/python /workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/scripts/campaign.py pilot

#!/usr/bin/env bash
set -euo pipefail
# Copy records only; never delete server or local historical records.
target_root=/home/lzy/projects/brava_flow_roi_18mlmin/population_v1
remote_root=/workspace/brava_flow_roi_18mlmin_20260928/population_v1
mkdir -p "$target_root/runs" "$target_root/analysis" "$target_root/validation/server" "$target_root/logs/server"
rsync -az "vast4090:$remote_root/runs/" "$target_root/runs/"
rsync -az "vast4090:$remote_root/analysis/" "$target_root/analysis/"
rsync -az "vast4090:$remote_root/validation/" "$target_root/validation/server/"
rsync -az "vast4090:$remote_root/logs/" "$target_root/logs/server/"
for pose in 0 15; do
  rsync -az "vast4090:$remote_root/candidate_$pose/" "$target_root/candidate_$pose/"
done
rsync -az "vast4090:$remote_root/figures/" "$target_root/figures/"

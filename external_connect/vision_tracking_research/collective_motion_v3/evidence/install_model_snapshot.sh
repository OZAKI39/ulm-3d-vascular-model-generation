#!/bin/bash
set -e
mkdir -p /workspace/f10w_vision_20261017/snapshots/77a2bf7e9
tar -xf /workspace/f10w_vision_20261017/setup/model_snapshot.tar -C /workspace/f10w_vision_20261017/snapshots/77a2bf7e9
cp /workspace/f10w_vision_20261017/setup/EXECUTION_models.json /workspace/f10w_vision_20261017/snapshots/77a2bf7e9/external_connect/vision_tracking_research/collective_motion_v3/EXECUTION.json
supervisorctl reread
supervisorctl update f10w_models_smoke_20261017 f10w_data_20261017
supervisorctl start f10w_models_smoke_20261017 f10w_data_20261017

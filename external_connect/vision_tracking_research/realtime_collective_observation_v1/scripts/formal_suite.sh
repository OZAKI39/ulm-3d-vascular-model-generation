#!/bin/bash
set -euo pipefail
export QT_QPA_PLATFORM=offscreen
run_root=/workspace/f10y_vision_20261017/runs/formal_v6
source_video=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete
mkdir -p "$run_root"
# Finish only our already-running stability pilot before timing this suite.
while supervisorctl status f10y_ui_stability_v3 | grep -q RUNNING; do sleep 5; done
python -m pip freeze > "$run_root/environment_remote.lock.txt"
python scripts/full_study.py --backend B4 --video "$source_video" --output "$run_root/study"
python scripts/synthetic_validation.py --backend B4 --video "$source_video" --output "$run_root/synthetic_frozen" --split frozen_test
python run_replay.py --backend B4 --video "$source_video" --autostart --seconds 33 --simulate-events --output "$run_root/ui" --capture-ui "$run_root/Video_F10Y_1_realtime_replay_observation.mp4"
python run_replay.py --backend B4 --video "$source_video" --autostart --loop --seconds 605 --output "$run_root/stability"

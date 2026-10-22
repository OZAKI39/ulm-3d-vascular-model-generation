#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export QT_QPA_PLATFORM=xcb
export PYOPENGL_PLATFORM=glx
export PYTHONDONTWRITEBYTECODE=1
mkdir -p logs
annotation_log="logs/annotation_$(date +%Y%m%dT%H%M%S)_$$.log"
printf '标注工作台正在启动；无需等待 SAM。日志：%s\n' "$PWD/$annotation_log"
printf '默认开启后台自动微调：仅使用人工确认的训练标签；界面可暂停。\n'
.venv_annotation/bin/python -u -B scripts/annotate.py "$@" 2>&1 | tee "$annotation_log"

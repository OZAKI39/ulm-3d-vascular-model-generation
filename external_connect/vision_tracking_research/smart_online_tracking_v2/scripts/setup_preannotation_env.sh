#!/usr/bin/env bash
# Optional Cellpose-SAM comparison environment; existing GUI / SAM env is unchanged.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
if [[ ! -x .venv_annotation/bin/python ]]; then
  printf '%s\n' 'Prepare the existing .venv_annotation environment first.' >&2
  exit 1
fi
if [[ ! -x .venv_preannotation/bin/python ]]; then
  .venv_annotation/bin/python -m venv .venv_preannotation
fi
.venv_preannotation/bin/python - <<'PY'
from pathlib import Path
import sysconfig,subprocess
root=Path.cwd()
shared=subprocess.check_output([str(root/'.venv_annotation/bin/python'),'-c','import sysconfig;print(sysconfig.get_paths()["purelib"])'],text=True).strip()
Path(sysconfig.get_paths()['purelib'],'annotation_dependencies.pth').write_text(shared+'\n')
PY
# Reuse the pinned GUI environment's torch/scipy/numpy without copying GPU libraries.
.venv_preannotation/bin/python -m pip install --no-deps \
  cellpose==4.2.1 fill-voids==2.1.2 roifile==2026.9.22 \
  opencv-python-headless==5.0.0.93 segment_anything==1.0

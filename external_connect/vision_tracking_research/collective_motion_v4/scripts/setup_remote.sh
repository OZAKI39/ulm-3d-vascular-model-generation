#!/bin/bash
set -euo pipefail
base=/workspace/f10x_piv_20261017
mkdir -p "$base"/{setup,logs,runs,third_party,snapshots}
test ! -d "$base/.venv"
/workspace/f10w_vision_20261017/.venv/bin/python -m venv "$base/.venv"
cp -a /workspace/f10w_vision_20261017/.venv/lib/python3.12/site-packages/. "$base/.venv/lib/python3.12/site-packages/"
source "$base/.venv/bin/activate"
python -m pip install --upgrade 'openpiv==0.26.1' maturin
git clone --depth 1 https://github.com/OpenPIV/openpiv-python.git "$base/third_party/openpiv-python"
git clone --depth 1 https://github.com/NikNazarov/TorchPIV.git "$base/third_party/TorchPIV"
git clone --depth 1 https://github.com/ali-sh-96/openpiv-python-gpu.git "$base/third_party/openpiv-python-gpu"
python -m pip freeze > "$base/setup/requirements_initial.lock"
python - <<'PY'
import importlib.metadata as m,inspect,json,torch,cupy,openpiv
from openpiv import pyprocess
report={'openpiv':m.version('openpiv'),'torch':torch.__version__,'cuda':torch.version.cuda,'cupy':cupy.__version__,'gpu':torch.cuda.get_device_name(),'signature':str(inspect.signature(pyprocess.extended_search_area_piv))}
try:
 import openpiv_rust
 report['rust_extension']=openpiv_rust.__file__
except Exception as e: report['rust_error']=repr(e)
print(json.dumps(report,indent=2))
open('/workspace/f10x_piv_20261017/setup/import_audit.json','w').write(json.dumps(report,indent=2))
PY

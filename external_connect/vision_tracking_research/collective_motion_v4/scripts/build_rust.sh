#!/bin/bash
set -euo pipefail
base=/workspace/f10x_piv_20261017
source "$base/.venv/bin/activate"
export CARGO_HOME="$base/toolchains/cargo" RUSTUP_HOME="$base/toolchains/rustup"
export PATH="$CARGO_HOME/bin:$PATH"
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs -o "$base/setup/rustup-init.sh"
sh "$base/setup/rustup-init.sh" -y --no-modify-path --profile minimal
cd "$base/third_party/openpiv-python"
CARGO_BUILD_JOBS=4 maturin develop --release -m crates/openpiv_rust/Cargo.toml
python -m pip install --no-deps -e .
python -m pip install --upgrade 'numpy==2.2.6' 'pandas>=2.2' 'scikit-learn>=1.5' 'pyarrow==19.0.1'
python - <<'PY'
import numpy as np,openpiv_rust,openpiv,torch,cupy,inspect,json,importlib.metadata as m
from openpiv import pyprocess
import cv2,av,pandas,sklearn,scipy,pyarrow
a=np.random.default_rng(10).normal(size=(192,192)).astype('float32')
u,v,s=pyprocess.extended_search_area_piv(a,np.roll(a,2,axis=1),48,overlap=24,backend='rust',normalized_correlation=True)
r={'openpiv':m.version('openpiv'),'rust_path':openpiv_rust.__file__,'signature':str(inspect.signature(pyprocess.extended_search_area_piv)),'rust_actual_output_median_uv':[float(np.median(u)),float(np.median(v))],'torch':torch.__version__,'cupy':cupy.__version__,'numpy':np.__version__,'gpu':torch.cuda.get_device_name()}
print(json.dumps(r,indent=2));open('/workspace/f10x_piv_20261017/setup/rust_actual_execution.json','w').write(json.dumps(r,indent=2))
PY
python -m pip freeze > "$base/setup/requirements_backend.lock"

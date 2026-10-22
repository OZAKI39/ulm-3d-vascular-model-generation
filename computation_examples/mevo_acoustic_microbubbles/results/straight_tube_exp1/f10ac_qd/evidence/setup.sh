#!/bin/bash
set -euo pipefail
cd /workspace/mevo_acoustic_microbubbles/f10ac_ad_20261018
test -x .venv/bin/python || uv venv --python /usr/bin/python3.12 .venv
uv pip install --python .venv/bin/python numpy==2.2.6 scipy==1.15.3 h5py==3.16.0 /workspace/mevo_acoustic_microbubbles/migration/wheels/warp_lang-1.18.0+cu12-py3-none-manylinux_2_28_x86_64.whl osqp==1.0.5 ribs==0.12.0 pymoo matplotlib==3.11.2 pillow==13.0.0 av==19.0.1 pymupdf==1.28.2 pyyaml pytest
uv pip freeze --python .venv/bin/python > logs/requirements.txt
git clone --quiet https://github.com/polycfd/apecss.git third_party/apecss
git -C third_party/apecss checkout a0ee071ca247b5489e997d90e83b5f53c40e4ec7
export APECSS_DIR=$PWD/third_party/apecss USRLIB_DIR=/usr/lib/x86_64-linux-gnu
cmake -S third_party/apecss/lib -B third_party/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-march=native
cmake --build third_party/build --parallel 2
gcc -O3 -march=native -DNDEBUG -I third_party/apecss/include third_party/apecss/examples/ultrasound/src/ultrasound_apecss.c third_party/build/libapecss.a -lm -o third_party/build/ultrasound_apecss
sha256sum third_party/build/ultrasound_apecss > logs/apecss.sha256

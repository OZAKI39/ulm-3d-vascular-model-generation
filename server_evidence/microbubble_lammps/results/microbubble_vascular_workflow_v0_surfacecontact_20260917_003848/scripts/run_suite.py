#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json
S=Path(__file__).resolve().parents[1]
assert 'surfacecontact_' in str(S) and str(S).startswith('/workspace/')
assert json.loads((S/'validation/SYNTHETIC_SURFACE_QUALIFICATION.json').read_text())['status']=='PASS'
assert json.loads((S/'validation/PRODUCTION_SURFACE_CPP_PYTHON.json').read_text())['status']=='PASS'
for name in ['LOW','MEDIUM','LONG_TRANSPORT']:
    for ranks in [1,4]:
        subprocess.run([sys.executable,str(S/'scripts/run_case.py'),name,str(ranks)],check=True)

from pathlib import Path
import csv, hashlib, json, sys
import numpy as np

HERE=Path(__file__).resolve().parent
FLOW=HERE.parent
REPORT=FLOW/'reports/taylor_hood_p2p1_validation_v1'
CTX=json.loads((REPORT/'run_context.json').read_text())
BASE=Path(CTX['baseline']); CASE=Path(CTX['case']); SOURCE=Path(CTX['source'])
PARTICLE=Path('/home/lzy/projects/ulm_particle_3d_particle0/particle_3d')
OLD=PARTICLE/'reports/particle9a3_interior_inlet'
OLD_DIAG=PARTICLE/'reports/particle9a3b_flowfield_conservation'
sys.path.insert(0,str(PARTICLE/'src'))

def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def csvwrite(path,rows):
    rows=list(rows)
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def baseline_arrays():
    a=np.load(BASE/'frozen_flow/flow_arrays_si.npz')
    return {k:a[k] for k in a.files}

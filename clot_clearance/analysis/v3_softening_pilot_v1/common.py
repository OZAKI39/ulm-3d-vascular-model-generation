from pathlib import Path
import json,hashlib,os,numpy as np
C=json.loads(Path(__file__).with_name('config.json').read_text());TASK=C['task'];RUN=C['run']
R=Path('/workspace/y_clot/project');A=R/'analysis/v3_softening_pilot_v1';O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
BASELINE=R/'verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T223001Z'
SMALL=BASELINE.parent/'20261014T210608Z';NATIVE=R/'verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z'
OLDV3=R/'verification/lps_final_native_bugfix_full_v3_regression_v1/20261011T211309Z'
HISTORY=R/'verification/y_clot_point_cloud_and_early_flow_coupling_v1/cases/C1_GAUSSIAN_A200'
OLDLIB=R/'verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z/build_owner_fix/liblammps.so'
LIB=O/'build/liblammps.so'
DX=50e-6;V=DX**3;RHO=1000.;MASS=RHO*V;G=940.;K=1e5/3;HORIZON=150.5e-6
BASE=.25e-6;MIN_H=BASE/8;END=.00233;CHECKPOINT_TRIGGER=.075
S0=np.array([0,.15,.140625,1.]);SC=np.array([0,.18,.16875,1.2])
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for q in iter(lambda:f.read(4194304),b''):h.update(q)
 return h.hexdigest()
def cv(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,Path):return str(x)
 raise TypeError(type(x).__name__)
def dump(p,z):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp')
 with tmp.open('w') as f:json.dump(z,f,indent=2,ensure_ascii=False,default=cv,allow_nan=False);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def save_state(p,s,**extra):
 p=Path(p);tmp=p.with_name(p.name+'.tmp')
 with tmp.open('wb') as f:np.savez_compressed(f,**s,**extra);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)

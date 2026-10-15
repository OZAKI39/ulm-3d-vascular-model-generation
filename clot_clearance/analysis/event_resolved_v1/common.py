from pathlib import Path
import hashlib,json,numpy as np,os
R=Path('/workspace/y_clot/project');TASK='lps_native_softening_event_resolved_integration_v1';RUN='20261013T235945Z'
A=R/'analysis/event_resolved_v1';O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
OLD=R/'verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z'
PREV=R/'verification/lps_native_softening_full_protocol_timestep_admission_v1/20261013T224838Z'
LIB=O/'build/liblammps.so';DX=50e-6;V=DX**3;RHO=1000.;G=940.;K=1e5/3;HORIZON=150.5e-6;S0=.15;SC=.18
BASE=0.25e-6;MIN_H=BASE/8;CHECKPOINT_TRIGGER=.075
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def dump(p,z):
 def cv(x):
  if isinstance(x,np.ndarray):return x.tolist()
  if isinstance(x,np.generic):return x.item()
  if isinstance(x,Path):return str(x)
  raise TypeError(type(x).__name__)
 Path(p).write_text(json.dumps(z,indent=2,ensure_ascii=False,default=cv,allow_nan=False))

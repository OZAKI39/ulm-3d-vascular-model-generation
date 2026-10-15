from pathlib import Path
import hashlib,json,numpy as np,os
R=Path('/workspace/y_clot/project');TASK='lps_native_softening_domain_retention_and_failure_safe_output_v1';RUN='20261014T150825Z'
A=R/'analysis/domain_retention_v1';O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
PREV=R/'verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z'
LIB=PREV/'build/liblammps.so';DX=50e-6;V=DX**3;RHO=1000.;G=940.;K=1e5/3;HORIZON=150.5e-6;S0=.15;SC=.18
BASE=.25e-6;MIN_H=BASE/8;CHECKPOINT_TRIGGER=.075;BOX=400e-6
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
 target=Path(p);temp=target.with_name(target.name+'.tmp')
 with temp.open('w') as f:json.dump(z,f,indent=2,ensure_ascii=False,default=cv,allow_nan=False);f.flush();os.fsync(f.fileno())
 os.replace(temp,target)
def save_state(path,state,**extra):
 path=Path(path);temp=path.with_name(path.name+'.tmp')
 with temp.open('wb') as f:np.savez_compressed(f,**state,**extra);f.flush();os.fsync(f.fileno())
 os.replace(temp,path)

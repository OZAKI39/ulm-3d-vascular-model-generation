"""Read-only baseline access and shared SI conventions for the mechanism audit."""
from pathlib import Path
import os,json,hashlib,importlib.util,numpy as np
R=Path(os.environ.get('Y_CLOT_PROJECT_ROOT','/workspace/y_clot/project'))
TASK='v3_chunk_fragmentation_mechanism_v1';RUN='20261015T174630Z'
A=R/'analysis'/TASK;O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
B=R/'verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z'
DX=50e-6;V=DX**3;RHO=1000.;G=940.;K=1e5/3;HORIZON=150.5e-6
S0=np.array([0,.15,.140625,1.]);SC=np.array([0,.18,.16875,1.2])
LIB=B/'build/liblammps.so'
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def cv(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,Path):return str(x)
 raise TypeError(type(x).__name__)
def dump(path,data):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.tmp')
 tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=cv,allow_nan=False)+'\n');os.replace(tmp,path)
def load_npz(path):
 with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def reference():return load_npz(B/'provenance/reference.npz')
def energy_model(ref):
 spec=importlib.util.spec_from_file_location('frozen_endpoint_star_model',R/'analysis/softening_feasibility_v1/model.py')
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod.SofteningLPS(ref,K=K,G=G)
def keyframes():
 for p in sorted((B/'full/keyframes').glob('step_*.npz')):yield p,load_npz(p)

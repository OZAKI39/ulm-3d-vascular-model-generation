from pathlib import Path
import os,json,hashlib,importlib.util
import numpy as np
R=Path(os.environ.get('Y_CLOT_PROJECT_ROOT',str(Path(__file__).resolve().parents[3])))
TASK='clot_irregular_lps_discretization_v4';RUN='20261019T221000Z'
A=R/'analysis'/TASK/'code';D=R/'analysis'/TASK;O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
V3_A=R/'analysis/clot_fragment_surface_coupling_v3';V3_O=R/'verification/clot_fragment_surface_coupling_v3/20261018T170400Z'
V2_A=R/'analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2';V2_O=R/'verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z'
LIB=O/'build/liblammps.so'
DX=50e-6;V=DX**3;RHO=1000.;K=1e5/3;G=940.;HORIZON=150.5e-6
BASE=.25e-6;MIN_H=BASE/8;END=.0036;RISE=.002;HOLD=.0024;OFF=.0028;COHORT=.0034
MU=.001;H_EFF=25e-6;BETA=MU/H_EFF
S0=np.array([0,.15,.140625,1.]);SC=np.array([0,.18,.16875,1.2])
LABEL='MECHANICAL_SHEAR_SURROGATE';HET_LABEL='EXPLORATORY_UNCALIBRATED'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def cv(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,Path):return str(x)
 raise TypeError(type(x).__name__)
def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);temp=p.with_suffix(p.suffix+'.tmp');temp.write_text(json.dumps(x,indent=2,default=cv,allow_nan=False)+'\n');os.replace(temp,p)
def load_npz(p):
 with np.load(p,allow_pickle=False) as f:return dict(f)
def module(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def energy_model(ref):return module('v4_energy',A/'energy.py').SofteningLPS(ref,K=K,G=G)

"""Isolated multifragment breakup; accepted historical evidence is read-only."""
from pathlib import Path
import os,json,hashlib,importlib.util
import numpy as np
R=Path(os.environ.get('Y_CLOT_PROJECT_ROOT','/workspace/y_clot/project'))
TASK='clot_biphasic_heterogeneous_shear_fragmentation_v2';RUN='20261017T231900Z'
A=R/'analysis'/TASK;O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
HISTORY=R/'verification/clot_multifragment_crushing_mechanism_v1/20261016T222934Z'
HA=R/'analysis/clot_multifragment_crushing_mechanism_v1'
BRACKETING=R/'verification/clot_mechanical_load_requirement_bracketing_v1/20261016T202533Z'
B=R/'verification/clot_multifragment_breakup_v1/20261017T110849Z';OLDLIB=B/'build/liblammps.so';LIB=Path(os.environ.get('HET_LIB',str(O/'build/liblammps.so')))
INV=R/'verification/physics_constrained_load_inverse_design_v1/20261016T144922Z'
CFD=R/'verification/y_clot_present_cfd_baseline_v1/20261011T154854Z'
MAP=R/'verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z'
DX=50e-6; V=DX**3; RHO=1000.; G=940.; K=1e5/3; HORIZON=150.5e-6
S0=np.array([0,.15,.140625,1.]); SC=np.array([0,.18,.16875,1.2]); BASE=.25e-6
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
def dump(p,z):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
 tmp.write_text(json.dumps(z,ensure_ascii=False,indent=2,default=cv,allow_nan=False)+'\n');os.replace(tmp,p)
def load_npz(p):
 with np.load(p,allow_pickle=False) as z:return dict(z)
def module(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def energy_model(ref):return module('frozen_energy',R/'analysis/softening_feasibility_v1/model.py').SofteningLPS(ref,K=K,G=G)

MIN_H=BASE/8
END=.0036

PREV=R/'verification/clot_acoustic_streaming_shear_multifragment_v1/20261017T132706Z'
PA=R/'analysis/clot_acoustic_streaming_shear_multifragment_v1'
MU=0.001; H_EFF=25e-6; BETA=MU/H_EFF
RISE=.002; HOLD=.0024; OFF=.0028; COHORT=.0034
CASES=['uniform','biphasic']
LABEL='PRESCRIBED_SHEAR_MECHANICAL_SURROGATE'

OLD=R/'verification/clot_shear_load_inverse_design_v1/20261017T165129Z'
OLD_A=R/'analysis/clot_shear_load_inverse_design_v1'
W_TARGET=25.698197e-9

MV=R/'verification/clot_shear_multivortex_fragmentation_v2/20261017T175456Z'
MV_A=R/'analysis/clot_shear_multivortex_fragmentation_v2'
HET_LABEL='EXPLORATORY_UNCALIBRATED_BIPHASIC_MATERIAL'

PREVIOUS=R/'verification/clot_spatially_heterogeneous_variational_lps_v1/20261017T214917Z'
PREVIOUS_A=R/'analysis/clot_spatially_heterogeneous_variational_lps_v1'
OLDLIB=PREVIOUS/'build/liblammps.so'
SEED=20261018
COUNTS=(25,25,13)

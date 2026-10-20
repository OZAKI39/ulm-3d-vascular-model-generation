"""Isolated inverse-design paths; accepted native assets are read-only."""
from pathlib import Path
import os, json, hashlib, importlib.util
import numpy as np
R=Path(os.environ.get('Y_CLOT_PROJECT_ROOT','/workspace/y_clot/project'))
TASK='clot_mechanical_load_requirement_bracketing_v1'; RUN='20261016T202533Z'
A=R/'analysis'/TASK; O=R/'verification'/TASK/RUN; F=R/'results'/TASK/RUN
B=R/'verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z'
PREV=R/'verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z'
SPEC=PREV/'benchmark_corrected_v3/shear_x/provenance'
CFD=R/'verification/y_clot_present_cfd_baseline_v1/20261011T154854Z'
MAP=R/'verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z'
LIB=B/'build/liblammps.so'
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

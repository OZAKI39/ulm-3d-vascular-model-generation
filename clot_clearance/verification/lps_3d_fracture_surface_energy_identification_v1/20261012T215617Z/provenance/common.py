from pathlib import Path
import json,hashlib,time,os
import numpy as np
ROOT=Path(__file__).resolve().parents[2];RUN=os.environ.get('Y_CLOT_DUAL_RUN_ID','20261011T154854Z')
COMMON=ROOT/'verification/y_clot_cfd_and_wide_pvce_dual_track_v1'/RUN
A=ROOT/'verification/lps_wide_pure_shear_pvce_feasibility_v1'/RUN
B=ROOT/'verification/y_clot_present_cfd_baseline_v1'/RUN
FA=ROOT/'results/lps_wide_pure_shear_pvce_feasibility_v1'/RUN
FB=ROOT/'results/y_clot_present_cfd_baseline_v1'/RUN
OLD=ROOT/'verification/lps_pvce_energy_release_feasibility_v1/20261011T132103Z'
V2=ROOT/'verification/lps_pvce_continuum_and_near_onset_audit_v2/20261011T142003Z'
V3=ROOT/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1'
DX=50e-6;DELTA=150.5e-6;G=940.;K=100000/3;RHO=1000.
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,v):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 def cv(o):
  if isinstance(o,np.ndarray):return o.tolist()
  if isinstance(o,np.generic):return o.item()
  if isinstance(o,Path):return str(o)
  raise TypeError(type(o).__name__)
 Path(p).write_text(json.dumps(v,indent=2,default=cv,allow_nan=False))

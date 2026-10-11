"""New-run paths and read-only imports; production code is never changed."""
from analysis.lps_pvce_v1.core import ROOT,DX,DELTA,G,K,RHO,model
from pathlib import Path
import numpy as np,json,hashlib,time
RUN='20261011T142003Z'
BASE=ROOT/'verification/lps_pvce_continuum_and_near_onset_audit_v2'/RUN
OLD=ROOT/'verification/lps_pvce_energy_release_feasibility_v1/20261011T132103Z'
FIG=ROOT/'results/lps_pvce_continuum_and_near_onset_audit_v2'/RUN
NAMES=['NOTCH_INTERIOR_FIRST','NOTCH_SIDE_SURFACE_FIRST','GRIP_FIRST','OTHER_SURFACE_FIRST','FREE_INTERIOR_FIRST']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def regions(X,b,tip=-100e-6):
 m=X[b].mean(1);dt=np.hypot(m[:,0]-tip,m[:,2]);dg=np.minimum(abs(m[:,2]-350e-6),abs(m[:,2]+350e-6));ds=np.min(np.array([300,300,500])*1e-6-abs(m),axis=1)
 # Identical V2/V3 notch>grip>surface>interior priority; notch split by y-side distance.
 reg=np.full(len(b),4,dtype=np.int8);reg[ds<=DELTA]=3;reg[dg<=DELTA]=2;reg[dt<=DELTA]=0;reg[(dt<=DELTA)&(300e-6-abs(m[:,1])<=DELTA)]=1
 return reg

def risk(z):
 b=z['bonds'];X=z['X'];x=z['x_m'];alive=z['alive'].astype(bool);r=z['r'];reg=regions(X,b)
 st=(np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)-r)/r;ids=np.flatnonzero(alive);k=ids[np.argmax(st[ids])];mx=float(st[k]);rm={n:float(st[alive&(reg==j)].max()/.15) for j,n in enumerate(NAMES)}
 # Competing event if distinct spatial classes lie within 2% of global utilization.
 competitors=[n for n,v in rm.items() if v>=.98*mx/.15];label='COMPETING_FIRST' if len(competitors)>1 else NAMES[int(reg[k])]
 near=ids[st[ids]>=mx*(1-1e-5)]
 return dict(max_stretch=mx,max_utilization=mx/.15,location=label,highest_region=NAMES[int(reg[k])],competing_regions=competitors,regional_utilization=rm,maximum_reference_bond_id=int(k+1),maximum_atom_ids=(b[k]+1).tolist(),maximum_midpoint_um=(X[b[k]].mean(0)*1e6).tolist(),symmetric_tie_reference_ids=(near+1).tolist(),symmetric_tie_atom_ids=(b[near]+1).tolist(),symmetric_tie_midpoints_um=(X[b[near]].mean(1)*1e6).tolist())

def save(p,obj):
 def cv(v):
  if isinstance(v,np.ndarray):return v.tolist()
  if isinstance(v,np.generic):return v.item()
  raise TypeError(type(v).__name__)
 Path(p).write_text(json.dumps(obj,indent=2,default=cv,allow_nan=False))

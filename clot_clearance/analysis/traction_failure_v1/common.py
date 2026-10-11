"""Isolated offline audit. No production model, source field or topology mutation."""
from pathlib import Path
import hashlib,json,os,time,resource
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
C=json.loads((Path(__file__).parent/'config.json').read_text())
RUN=C['run_id'];COMMON=ROOT/C['common'];A=ROOT/C['A'];B=ROOT/C['B'];FA=ROOT/C['FA'];FB=ROOT/C['FB']
PREV=ROOT/'verification/y_clot_present_cfd_baseline_v1/20261011T154854Z'
V3=ROOT/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1'
M5C=ROOT/'verification/lammps_notch_fracture_grip_transfer_v3/20261011T115626Z/M5C'
DX=50e-6;DELTA=150.5e-6;V=DX**3;K=33333.333333333336;G=940.
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def convert(v):
 if isinstance(v,np.ndarray):return v.tolist()
 if isinstance(v,np.generic):return v.item()
 if isinstance(v,Path):return str(v)
 raise TypeError(type(v).__name__)
def save(p,v):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,default=convert,allow_nan=False)+'\n')
def stats(a):
 a=np.asarray(a);return dict(zip(['min','p05','median','p95','max'],np.quantile(a,[0,.05,.5,.95,1]))) if a.size else {'count':0}
def cost(t,c):return dict(wall_s=time.perf_counter()-t,process_cpu_s=time.process_time()-c,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,gpu_compute=False)
def closest_surface(points,triangles):
 import vtk,pyvista as pv
 poly=pv.PolyData(triangles.reshape(-1,3),np.c_[np.full(len(triangles),3),np.arange(len(triangles)*3).reshape(-1,3)])
 loc=vtk.vtkStaticCellLocator();loc.SetDataSet(poly);loc.BuildLocator();near=[];ids=[];dist=[]
 for point in points:
  cp=[0.,0.,0.];cid=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.);loc.FindClosestPoint(point,cp,cid,sub,d2);near.append(cp);ids.append(int(cid));dist.append(float(d2)**.5)
 return np.array(near),np.array(ids),np.array(dist)

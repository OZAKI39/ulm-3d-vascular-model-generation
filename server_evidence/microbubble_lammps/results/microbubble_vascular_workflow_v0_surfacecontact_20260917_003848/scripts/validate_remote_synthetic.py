from pathlib import Path
import csv,json,collections,numpy as np
from flow_geometry import Geometry
from independent_wall import segment_distance
S=Path(__file__).resolve().parents[1]
with (S/'raw/REMOTE_SYNTHETIC_TRAJECTORIES.csv').open() as f:rows=list(csv.DictReader(f))
groups=collections.defaultdict(list)
for r in rows:groups[(r['mesh'],float(r['dt_max']))].append(r)
local=json.loads((S/'raw/SYNTHETIC_TRAJECTORIES.json').read_text());out=[]
for (name,dt),rr in groups.items():
 d=np.load(S/'raw/synthetic'/f'{name}.npz');path=S/'raw/synthetic'/f'{name}_vtk.npz';np.savez(path,points=d['vertices'],faces=d['faces'],classes=np.zeros(len(d['faces']),int));g=Geometry(path);tr=np.array([[float(r[k]) for k in ['time','x','y','z','dt','stage']] for r in rr]);end=tr[tr[:,-1]==1]
 for r in rr:assert float(r['swept_gap'])>=1e-10-3e-20 and float(r['correction'])<=2.5e-11
 if name.startswith('plane'):
  expected_time=tr[:,0]-np.where(tr[:,-1]==0,tr[:,4]/2,0);assert np.max(abs(tr[:,1]-(-1e-6+1e-4*expected_time)))<1e-18;assert np.max(abs(tr[:,2]))<1e-18;assert np.max(abs(tr[:,3]-(1e-6+1e-10+1e-17)))<1e-18
 else:
  a=next(a for a in local if a['name']==name and a['dt']==dt and abs(a['time']-.01)<1e-12);ref=np.array(a['trajectory']);assert len(end)==len(ref);assert np.max(abs(end[:,0:4]-ref[:,0:4]))<1e-18
 # Independent continuous triangle check on every remote stage after the initial step.
 minimum=1.
 for k in range(2,len(tr),2):
  base=tr[k-1,1:4]
  for j in [k,k+1]:
   distance=segment_distance(g,base,tr[j,1:4],1e-6+1e-10);minimum=min(minimum,distance-1e-6);assert distance>=1e-6+1e-10-2e-14
 out.append(dict(mesh=name,dt=dt,stages=len(rr),status='PASS',minimum_independent_swept_gap=minimum,final_time=float(end[-1,0])))
(S/'validation/REMOTE_NATIVE_SYNTHETIC_INDEPENDENT.json').write_text(json.dumps(dict(status='PASS',cases=out,checks='Remote C++ trajectories vs local independently checked fixtures, analytic plane and full NumPy/VTK swept safety'),indent=2)+'\n');print('REMOTE_NATIVE_SYNTHETIC_INDEPENDENT_PASS',len(out))

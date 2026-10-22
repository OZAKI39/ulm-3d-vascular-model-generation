"""Separate global distribution shape from actual local geometric components.

Read-only postprocessing. Threshold sensitivity is diagnostic, never a change to
archive quality/admission. Initial shapes use each final component's fixed IDs.
"""
import time
import h5py
import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components
from scripts.straight_tube_exp1.f10ab.metrics import geometry
from .common import *

def partition(x0,x,factor):
 d=cdist(x,x);nc,labels=connected_components(d<=factor*R,directed=False)
 center=x.mean(0);groups=[];within=0.;between=0.
 for label in range(nc):
  ids=np.flatnonzero(labels==label);y=x[ids];c=y.mean(0);rg=float(np.sqrt(np.mean(np.sum((y-c)**2,axis=1))));weight=len(ids)/len(x)
  r0=float(np.sqrt(np.mean(np.sum((x0[ids]-x0[ids].mean(0))**2,axis=1))))
  item=dict(particle_ids=ids,size=len(ids),center_m=c,Rg_m=rg,initial_same_ID_Rg_m=r0,Rg_fraction_change=rg/r0-1 if r0>0 else None)
  if len(ids)>=4:
   g=geometry(y);item.update(aspect_3d=g['aspect_3d'],principal_rms_lengths_m=g['principal_rms_lengths_3d_m'],orientation_valid=g['orientation_valid_3d'])
  else:item.update(aspect_3d=None,shape_scope='Fewer than four particles; no full 3D aspect admission')
  groups.append(item);within+=weight*rg**2;between+=weight*np.sum((c-center)**2)
 total=float(np.mean(np.sum((x-center)**2,axis=1)));assert np.isclose(total,within+between,rtol=1e-12,atol=1e-30)
 groups.sort(key=lambda g:(-g['size'],int(g['particle_ids'][0])))
 distances=[float(d[np.ix_(a['particle_ids'],b['particle_ids'])].min()) for i,a in enumerate(groups) for b in groups[i+1:]]
 return dict(threshold_R=factor,components=nc,groups=groups,global_Rg_m=np.sqrt(total),within_component_Rg_squared_m2=within,between_component_Rg_squared_m2=between,minimum_intercomponent_center_distance_m=min(distances) if distances else None)

def main():
 start=time.perf_counter();rows=read(BASE/'ALL_RUNS.json');out=[]
 for r in rows:
  p=BASE/'dataset'/(r['id']+'.h5')
  with h5py.File(p) as h:x0=h['positions_m'][0];x=h['positions_m'][-1];assert h.attrs['raw_sha256']==r['trajectory_sha256']
  analyses=[partition(x0,x,f) for f in [2.05,2.1,2.15]]
  assert analyses[1]['components']==r['final']['proximity_components']
  out.append(dict(id=r['id'],scientific_status=r['scientific_status'],completed_s=r['completed_s'],raw_sha256=r['trajectory_sha256'],endpoint_analyses=analyses))
 write(BASE/'GROUP_DISTRIBUTION_DIAGNOSTICS.json',dict(rows=out,threshold_policy='2.1R remains the original admission rule; 2.05R and 2.15R diagnostic only',physical_interpretation='Proximity is not adhesion or fusion; group identities may change',wall_s=time.perf_counter()-start,identity=identity()))
 chosen=read(BASE/'ELITE_VALIDATION_PREREGISTRATION.json')['selected']+[read(BASE/'OUTSIDE_PREREGISTRATION.json')['chosen']];lookup={r['id']:r for r in out}
 lines=['# Whole distribution and local components','','The complete fixed-ID distribution and its largest geometric component are evaluated separately. Initial local Rg uses the same IDs as the final local component. Changing component membership is not physical fusion.','','| Case | Global Rg (um) | Components | Largest group N | Largest group Rg (um) | Same-ID local Rg change (%) | Components at 2.05/2.1/2.15R |','| --- | --- | --- | --- | --- | --- | --- |']
 for r in chosen:
  a=lookup[r['id']]['endpoint_analyses'];z=a[1];g=z['groups'][0];lines.append(f"| {r['id']} | {z['global_Rg_m']*1e6:.5f} | {z['components']} | {g['size']} | {g['Rg_m']*1e6:.5f} | {g['Rg_fraction_change']*100:.4f} | {' / '.join(str(d['components']) for d in a)} |")
 lines+=['','All 366 retained endpoints are in GROUP_DISTRIBUTION_DIAGNOSTICS.json. Incomplete-prefix endpoints are labeled and are not admitted full trajectories. Weighted within-component and between-component variances independently reconstruct total Rg squared.','', 'EXPLORATORY | PHYSICS_UNCALIBRATED | training_ready=false']
 (BASE/'LOCAL_COMPONENT_SHAPES.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()

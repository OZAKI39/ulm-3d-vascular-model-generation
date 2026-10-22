"""Real native P1 streamlines and local velocity glyphs for the BraVa case."""
from pathlib import Path
import sys,json,hashlib,os
from types import SimpleNamespace
from collections import Counter
import numpy as np,pyvista as pv
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path(os.environ.get('BRAVA_PARTICLE_ROOT','/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms'))
sys.path.insert(0,str(SOURCE/'particle_3d/src'))
from particle_3d.field import FrozenFEMField
from particle_3d.inlet_flux import InletFluxSampler
from particle_3d.validation_boundary import ValidationBoundaryClassifier
from particle_3d.particle82_point_native import NativePointTracer
active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text());case=ROOT/active['case'];out=ROOT/'visualization/data';out.mkdir(parents=True,exist_ok=True)
a=np.load(case/'frozen_flow/flow_arrays_si.npz');field=FrozenFEMField(a['points_m'],a['tetra'],a['velocity_m_s'],a['pressure_pa'])
roles=['INLET','OUTLET_01','OUTLET_02','OUTLET_03'];boundaries={r:pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f'{r}.vtp') for r in roles}
classifier=ValidationBoundaryClassifier({r:s for r,s in boundaries.items() if r.startswith('OUTLET_')})
env=SimpleNamespace(field=field,boundaries=boundaries,classifier=classifier);native=NativePointTracer(env)
s=boundaries['INLET'];idx=np.asarray(s['GlobalNodeID'],int)[s.faces.reshape(-1,4)[:,1:]]-1
xyz=a['points_m'][idx];n=-np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
flux=np.einsum('nki,ni->nk',a['velocity_m_s'][idx],n);sampler=InletFluxSampler(xyz,flux)
seeds,_=sampler.sample(np.random.default_rng(2026092821),192)
records=[];paths={};points=[];velocities=[];lines=[];offset=0;selected=[]
for i,seed in enumerate(seeds):
 row=native.trace(seed,step_m=2.5e-5,error=2.5e-9,horizon_m=.6);path=row.pop('path');paths[f'seed_{i:03d}']=path
 records.append(dict(seed_id=i,**row))
 if row['outlet'] and len(selected)<96:
  v=np.empty((len(path),3));cell=field.locate(path[0,1:])[0]
  for j,p in enumerate(path[:,1:]):
   position=np.ascontiguousarray(p);nextcell=native.lib.p82_sample(position.ctypes.data,cell,v[j].ctypes.data)
   if nextcell<0:
    assert j==len(path)-1 and row['outlet'], 'Interior streamline sample outside its source field'
    v[j]=v[j-1]
   else:cell=nextcell
  selected.append(dict(seed_id=i,outlet=row['outlet'],transit_s=float(path[-1,0])))
  points.append(path[:,1:]);velocities.append(v);lines.extend([len(path),*range(offset,offset+len(path))]);offset+=len(path)
 if (i+1)%32==0:print('STREAMLINES',i+1,dict(Counter(r['outlet'] for r in records)),flush=True)
np.savez_compressed(out/'candidate_streamlines.npz',**paths)
(out/'candidate_trace_records.json').write_text(json.dumps(records,indent=2)+'\n')
assert len(selected)==96 and {r['outlet'] for r in selected}==set(roles[1:])
refinement=[]
for role in roles[1:]:
 for row in [r for r in selected if r['outlet']==role][:3]:
  coarse=paths[f"seed_{row['seed_id']:03d}"];fine=native.trace(coarse[0,1:],step_m=1.25e-5,error=6.25e-10,horizon_m=.6)
  assert fine['outlet']==role
  def resample(path):
   arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path[:,1:],axis=0),axis=1))]
   return np.column_stack([np.interp(np.linspace(0,arc[-1],1001),arc,path[:,j]) for j in (1,2,3)])
  error=float(np.linalg.norm(resample(coarse)-resample(fine['path']),axis=1).max())
  assert error<1e-5,(row,error)
  refinement.append(dict(seed_id=row['seed_id'],outlet=role,maximum_path_difference_m=error))

mesh=pv.PolyData(np.vstack(points),lines=np.array(lines));mesh['Velocity_m_s']=np.vstack(velocities);mesh['Speed_mm_s']=np.linalg.norm(mesh['Velocity_m_s'],axis=1)*1000;mesh.cell_data['Line_id']=np.arange(len(selected));mesh.save(out/'streamlines_si.vtp')
np.savez_compressed(out/'candidate_streamlines.npz',**paths)
# Centerline forks are identified by parent/child topology, then transformed, never copied from the mouse model.
swc=np.loadtxt(ROOT/'inputs/fitted_centerline.swc');ids=swc[:,0].astype(int);children={i:np.count_nonzero(swc[:,6]==i) for i in ids};forks=np.array([row[2:5] for row in swc if children[int(row[0])]>1])
T=np.array(json.loads((ROOT/'inputs/print_transform.json').read_text())['transform_4x4']);g=json.loads((ROOT/'inputs/geometry.json').read_text());origin=np.array(g['origin_source_mm'])
forks_m=(forks@T[:3,:3].T+T[:3,3]-origin)*1e-3;focus=forks_m[0]
t=a['tetra'];c=a['points_m'][t].mean(1);u=a['velocity_m_s'][t].mean(1);speed=np.linalg.norm(u,axis=1)
mask=(np.linalg.norm((c-focus)/np.array([.006,.006,.006]),axis=1)<.9)&(speed>1e-5)
indices=np.flatnonzero(mask);candidate=c[indices];chosen=[int(np.argmin(np.linalg.norm(candidate-focus,axis=1)))];d=np.linalg.norm(candidate-candidate[chosen[0]],axis=1)
while len(chosen)<350 and d.max()>.0006:
 j=int(np.argmax(d));chosen.append(j);d=np.minimum(d,np.linalg.norm(candidate-candidate[j],axis=1))
indices=indices[chosen];cloud=pv.PolyData(c[indices]);cloud['Velocity_m_s']=u[indices];cloud['Direction']=u[indices]/speed[indices,None];cloud['Speed_mm_s']=speed[indices]*1000;cloud['Tetra_id']=indices;cloud.save(out/'velocity_samples_si.vtp')
report=dict(data_sha256={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['streamlines_si.vtp','velocity_samples_si.vtp']},refinement_comparisons=refinement,source_flow_sha256=active['flow_sha256'],candidate_count=len(records),candidate_outcomes=dict(Counter(str(r['outlet']) for r in records)),selected_streamlines=selected,
 streamline_semantics='P1 RK23 massless flow lines, not finite-size microbubble dynamics; line counts are not flow fractions',step_m=2.5e-5,local_error_m=2.5e-9,
 selection='First 96 naturally completed source-order seeds; no path steering or extension',branch_focus_m=focus.tolist(),bifurcations_m=forks_m.tolist(),glyph_count=len(indices),glyph_direction='P1 velocity averaged at actual tetrahedron centroid',terminal_streamline_color='If exact cap crossing is outside the locator, repeat last interior speed for endpoint color only; path endpoint is unchanged')
(out/'PREPARATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='selected_streamlines'},indent=2))

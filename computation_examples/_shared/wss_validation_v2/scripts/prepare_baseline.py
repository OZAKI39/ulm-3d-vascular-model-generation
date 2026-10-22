"""Copy immutable H0 evidence, then qualify its saved solution for mesh comparisons."""
from pathlib import Path
import shutil,json
import numpy as np
import pyvista as pv
from case_common import *
from flow_solver_support.wss_case import coordinate_identity
H=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1')
case=V/'stage3/vessel_baseline';assert not case.exists();case.mkdir()
for name in ['SV_MESH','frozen_flow']:shutil.copytree(H/name,case/name)
(case/'run').mkdir()
for name in ['solver.xml','solver.log']:
 shutil.copy2(H/'run'/name,case/'run'/name)
shutil.copy2(H/'policy.json',case/'policy.json')
a=np.load(case/'SV_MESH/mesh_arrays.npz');x=a['points_m'];t=a['tetra'];b=a['boundary_triangles'];tags=a['facet_tags'];wall=b[tags==1];own=wss.boundary_owners(t,wall);cent,area,normal=wss.wall_geometry(x,t,wall,own);vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6
policy=json.loads((case/'policy.json').read_text());records=[];previous=None

def l2(v):
 arr=v[t];return float(np.sqrt(np.dot(vol,(np.sum(arr*arr,axis=tuple(range(1,arr.ndim)))+np.sum(arr.sum(axis=1)**2,axis=-1))/20))) if v.ndim==2 else float(np.sqrt(np.dot(vol,(np.sum(arr*arr,axis=1)+arr.sum(axis=1)**2)/20)))
for step in [10,20,30,40,50,60,70,71]:
 path=H/('run/1-procs/result_%03d.vtu'%step);g=pv.read(path);coordinate_identity(g.points,x);u=np.asarray(g['Velocity']);p=np.asarray(g['Pressure']).reshape(-1);ww=np.linalg.norm(wss.tangential_traction(wss.p1_gradients(x,t[own],u),normal,MU),axis=1);row=dict(step=step,source=str(path),sha256=sha(path),time_s=step*policy['dt_s'])
 if previous is not None:
  ou,op,ow=previous;row.update(velocity_relative_change=l2(u-ou)/l2(u),pressure_relative_change=l2(p-op)/max(l2(p),np.sqrt(vol.sum())),wss_area_L2_relative_change=float(np.sqrt(np.average((ww-ow)**2,weights=area)/np.average(ww**2,weights=area))))
  row['steady_interval_pass']=max(row[k] for k in ['velocity_relative_change','pressure_relative_change','wss_area_L2_relative_change'])<1e-7
 previous=(u.copy(),p.copy(),ww.copy());records.append(row)
assert all(r['steady_interval_pass'] for r in records[-4:])
# Original material, BC, mesh, and frozen solution are byte-identical copies.
paths=[p for root in [case/'SV_MESH',case/'frozen_flow',case/'run'] for p in root.rglob('*') if p.is_file()]+[case/'policy.json']
for p in paths:assert sha(p)==sha(H/p.relative_to(case))
dump(case/'reports/baseline_reuse.json',dict(status='QUALIFIED_REUSED_CFD',source_case=str(H),case_inputs_and_result_byte_identical=True,records=records,criteria='last 3 saved intervals u,p,WSS relative change <1e-7, separately check final mass and solver logs',note='No new CFD run; original independently solved H0 step71. Saved intervals are 10 steps except final 1 step.',source_hashes={str(p.relative_to(case)):sha(p) for p in paths}))
shutil.copy2(C/'reports/mesh_and_flow/geometry_qc.json',case/'reports/formal_geometry_qc.json')
shutil.copy2(C/'reports/mesh_and_flow/mesh_validity.json',case/'reports/formal_mesh_validity.json')
print(json.dumps(records[-4:],indent=2))

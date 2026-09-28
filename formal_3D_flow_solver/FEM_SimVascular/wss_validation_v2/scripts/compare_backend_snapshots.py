"""Compare explicitly identified execution configurations at a matched CFD step."""
import argparse,json
import numpy as np
import pyvista as pv
from case_common import *
from flow_solver_support.wss_case import coordinate_identity
p=argparse.ArgumentParser();p.add_argument('--reference','--gpu',dest='gpu',type=Path,required=True);p.add_argument('--candidate','--cpu',dest='cpu',type=Path,required=True);p.add_argument('--change-description',default='PETSc array backend and MPI partition differ');p.add_argument('--step',type=int,required=True);p.add_argument('--output',type=Path);a=p.parse_args();mg=np.load(a.gpu/'SV_MESH/mesh_arrays.npz');mc=np.load(a.cpu/'SV_MESH/mesh_arrays.npz')
for k in ['points_m','tetra','boundary_triangles','facet_tags']:assert np.array_equal(mg[k],mc[k])
x=mg['points_m'];t=mg['tetra'];b=mg['boundary_triangles'];wall=b[mg['facet_tags']==1];own=wss.boundary_owners(t,wall);_,area,n=wss.wall_geometry(x,t,wall,own);vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6
paths=[case/('run/%d-procs/result_%03d.vtu'%(json.loads((case/'policy.json').read_text()).get('MPI_ranks',1),a.step)) for case in [a.gpu,a.cpu]];grids=[pv.read(path) for path in paths]
for g in grids:coordinate_identity(g.points,x);assert np.array_equal(np.sort(g.cells_dict[pv.CellType.TETRA],axis=1),np.sort(t,axis=1))
def norm(z):
 q=z[t]
 if z.ndim==2:s=np.sum(q*q,axis=(1,2))+np.sum(q.sum(axis=1)**2,axis=1)
 else:s=np.sum(q*q,axis=1)+q.sum(axis=1)**2
 return float(np.sqrt(np.dot(vol,s/20)))
u,v=[np.asarray(g['Velocity']) for g in grids];p,q=[np.asarray(g['Pressure']).reshape(-1) for g in grids];wa,wb=[np.linalg.norm(wss.tangential_traction(wss.p1_gradients(x,t[own],uv),n,MU),axis=1) for uv in [u,v]]
r=dict(step=a.step,reference_case=a.gpu.name,candidate_case=a.cpu.name,source_snapshot_sha256=[sha(path) for path in paths],velocity_relative_L2=norm(v-u)/norm(u),pressure_relative_L2=norm(q-p)/norm(p),wss_area_relative_L2=float(np.sqrt(np.average((wb-wa)**2,weights=area)/np.average(wa**2,weights=area))),wss_max_abs_difference_Pa=float(np.max(abs(wb-wa))),comparison='same mesh, PDE, BC, timestep and iteration index; '+a.change_description+'; same linear residual targets',registered_relative_L2_limit=1e-9)
r['backend_equivalence_work_gate']=max(r[k] for k in ['velocity_relative_L2','pressure_relative_L2','wss_area_relative_L2'])<=1e-9;dump(a.output if a.output else V/'data'/('fine_backend_step%d.json'%a.step),r);print(json.dumps(r,indent=2))

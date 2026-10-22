"""Quantify the actual fine CFD correction from its recorded initial guess."""
from pathlib import Path
import argparse,json,hashlib,sys
import numpy as np
import pyvista as pv
from case_common import wss
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);a=p.parse_args();case=a.case.resolve();assert json.loads((case/'reports/execution.json').read_text())['status']=='PASS'
m=np.load(case/'SV_MESH/mesh_arrays.npz');f=np.load(case/'frozen_flow/flow_arrays_si.npz');initial=pv.read(case/'initial_state/flow_guess.vtu');x=m['points_m'];t=m['tetra'];b=m['boundary_triangles'];tags=m['facet_tags'];assert np.array_equal(initial.points,x)
u0=np.asarray(initial['Velocity']);p0=np.asarray(initial['Pressure']).ravel();u1=f['velocity_m_s'];p1=f['pressure_pa'];v=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6
wall=b[tags==1];own=wss.boundary_owners(t,wall);_,area,n=wss.wall_geometry(x,t,wall,own);mu=json.loads((case/'policy.json').read_text())['mu_Pa_s']
w0,w1=[np.linalg.norm(wss.tangential_traction(wss.p1_gradients(x,t[own],u),n,mu),axis=1) for u in [u0,u1]]
def norm(z):
 z=z[t];q=np.sum(z*z,axis=(1,2))+np.sum(z.sum(axis=1)**2,axis=1) if z.ndim==3 else np.sum(z*z,axis=1)+z.sum(axis=1)**2
 return float(np.sqrt(np.dot(v,q/20)))
ex=json.loads((case/'reports/execution.json').read_text());r=dict(case=case.name,actual_CFD_final_step=ex['final_step'],accepted_final_execution=True,initial_guess_is_not_final_solution=True,velocity_volume_L2_change_relative_to_final=norm(u1-u0)/norm(u1),pressure_volume_L2_change_relative_to_final=norm(p1-p0)/norm(p1),WSS_area_L2_change_relative_to_final=float(np.sqrt(np.average((w1-w0)**2,weights=area)/np.average(w1*w1,weights=area))),WSS_max_abs_change_Pa=float(np.max(abs(w1-w0))),initial_WSS_area_mean_Pa=float(np.average(w0,weights=area)),final_WSS_area_mean_Pa=float(np.average(w1,weights=area)),initial_guess_sha256=hashlib.sha256((case/'initial_state/flow_guess.vtu').read_bytes()).hexdigest(),final_field_sha256=hashlib.sha256((case/'frozen_flow/flow_arrays_si.npz').read_bytes()).hexdigest(),interpretation='Initial guess statistics are diagnostic only and excluded from mesh-convergence tables')
(case/'reports/initial_to_final_comparison.json').write_text(json.dumps(r,indent=2)+'\n');(case.parents[1]/'data/fine_initial_to_final_comparison.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

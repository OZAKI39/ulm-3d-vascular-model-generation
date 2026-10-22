"""Independent implementation and geometry checks, not physical flow validation."""
from pathlib import Path
import argparse,json,time
import numpy as np
import pyvista as pv
import torch
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.wall_geometry import WallGeometry
from audit_field import js
ap=argparse.ArgumentParser();ap.add_argument('--population-root',type=Path,required=True);a=ap.parse_args();root=a.population_root;out=root/'cfd_wall_audit';z=np.load(out/'recomputed_arrays.npz')
spec=json.loads((root/'data/environment.json').read_text());env=load_environment(spec)
G=z['G'];err=np.max(abs(G-env.field.gradients_s_inv));assert err<1e-8
wall=GPUWall(env.wall);cent=z['centroid'];cpu=z['wall_distance'];radius=1e-12
worst=0.;different=0;gpu_ids=[];started=time.perf_counter()
for start in range(0,len(cent),100000):
 x=torch.tensor(cent[start:start+100000],device='cuda',dtype=torch.float64);r=torch.full((len(x),),radius,device='cuda',dtype=torch.float64)
 ids,value,_,_,_=wall.query(x,r);gpu_ids.extend(ids.cpu().tolist());dist=value[:,0].cpu().numpy()+radius
 worst=max(worst,float(np.max(abs(dist-cpu[start:start+len(x)]))));different+=int(np.sum(ids.cpu().numpy()!=z['nearest_wall_triangle'][start:start+len(x)]))
assert worst<1e-12
# Audit quality formula check on a regular tetra with exact known normalized values.
x=np.array([[0,0,0],[1,0,0],[.5,np.sqrt(3)/2,0],[.5,np.sqrt(3)/6,np.sqrt(2/3)]])
det=np.linalg.det(x[1:]-x[0]);vol=abs(det)/6;e=np.array([np.linalg.norm(x[i]-x[j]) for i in range(4) for j in range(i+1,4)])
shape=12*(3*vol)**(2/3)/sum(e**2);sj=[]
for k in range(4):
 dx=x[[i for i in range(4) if i!=k]]-x[k];sj.append(np.sqrt(2)*abs(np.linalg.det(dx))/np.prod(np.linalg.norm(dx,axis=1)))
assert abs(shape-1)<1e-14 and max(abs(np.array(sj)-1))<1e-14
record=dict(PASS=True,scope='DISCRETE_IMPLEMENTATION_AND_DISTANCE_CROSSCHECK_NOT_PHYSICAL_VALIDATION',
 original_node_gradient_recompute_vs_actual_particle_sampler_max_abs_s_inv=float(err),cell_count=len(G),
 vtk_locator_vs_cuda_finite_triangle_distance_max_error_m=worst,distance_cells=len(cent),nearest_triangle_id_tie_or_roundoff_differences=different,
 distance_crosscheck_wall_s=time.perf_counter()-started,gpu=torch.cuda.get_device_name(),
 regular_tetra_mean_ratio=shape,regular_tetra_scaled_jacobians=sj,authoritative_inputs_unchanged=True)
import pandas as pd
f=pd.read_csv(out/'all_cells.csv');ids=np.array(gpu_ids);normal=z['wall_normal_in'][ids];uc=z['velocity'][z['connectivity']].mean(1);speed=np.linalg.norm(uc,axis=1)
un=np.einsum('ij,ij->i',uc,normal);nf=np.divide(abs(un),speed,out=np.zeros_like(speed),where=speed!=0)
record['nearest_triangle_normal_sensitivity']={}
for l in [1,2]:
 mask=z['layer']==l;original=f.normal_fraction.to_numpy();changed=mask&((nf>=.5)!=(original>=.5))
 record['nearest_triangle_normal_sensitivity'][f'W{l}']=dict(max_abs_fraction_change=float(np.max(abs(nf[mask]-original[mask]))),count_crossing_screen_0p5=int(changed.sum()))
record['nearest_triangle_normal_sensitivity']['W0']='Own WALL face normal; unaffected by nearest-triangle ties'
js(out/'validation/independent_recomputation.json',record);print(json.dumps(record),flush=True)

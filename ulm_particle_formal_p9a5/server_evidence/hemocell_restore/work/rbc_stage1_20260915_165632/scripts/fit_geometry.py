#!/usr/bin/env python3
"""Pre-solver deterministic placement audit using native HemoCell RBC vertices and frozen closed lumen."""
from pathlib import Path
import json,hashlib,time,csv,math,sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy,numpy_to_vtk,numpy_to_vtkIdTypeArray
C=Path(__file__).resolve().parent;G=C/'geometry';G.mkdir(exist_ok=True)
B=Path('/home/lzy/projects/compre_output/step3b_remote_bundle/20260913_122842/remote_bundle')
policy=dict(schema='PRE_SOLVER_GEOMETRY_SEARCH_V1',cell_scale_formula='cbrt(50/90)',cell_radius_m='3.91e-6*cbrt(50/90)',shape='native RBC_FROM_SPHERE',minimum_triangles=600,coarse_center_lattice_stride=2,normal_hemisphere_count=128,roll_degrees=[0,30,60],refinement_candidates=48,refinement_rounds=8,refinement_translation_initial_lu=1.,refinement_rotation_initial_deg=5.,CASE_A_C_min_wall_clearance_lu=2.,CASE_B_min_wall_clearance_lu=.5,penetration_tolerance_lu=0.,positive_signed_distance_means_outside=True,triangle_contact_test='vtkCollisionDetectionFilter, cell and box tolerance0, all contacts',geometry_metric_scope='Original closed STL wall plus caps; conservative relative to native outward inflate0.001LU',search_scope='Bounded deterministic rigid placement search, not a proof over all continuous poses; no mesh resizing, remeshing or timestep execution')
policyfile=C/'GEOMETRY_SEARCH_POLICY.json'
if policyfile.exists():assert json.loads(policyfile.read_text())==policy
else:policyfile.write_text(json.dumps(policy,indent=2)+'\n')
if '--freeze-only' in sys.argv:print('SEARCH_POLICY_FROZEN');sys.exit(0)
t0=time.monotonic();meta=json.loads((G/'NATIVE_MESH_GENERATION.json').read_text());dx=meta['dx_m'];dxum=dx*1e6
pg=json.loads((B/'frozen_inputs/step2/diagnostics/palabos_geometry.json').read_text());nx,ny,nz=pg['lattice_shape'];origin=np.array(pg['physical_origin_m']);stl=B/'frozen_inputs/step1_geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl'
verts=np.genfromtxt(G/'RBC_REFERENCE_VERTICES.csv',delimiter=',',skip_header=1)[:,1:];faces=np.genfromtxt(G/'RBC_REFERENCE_TRIANGLES.csv',delimiter=',',skip_header=1,dtype=int)[:,1:]
assert len(faces)>=600 and np.all(np.isfinite(verts));verts=verts-verts.mean(axis=0)
tri=verts[faces];vol=abs(np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6)*dxum**3;area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2*dxum**2
reader=vtk.vtkSTLReader();reader.SetFileName(str(stl));reader.Update();surface=vtk.vtkPolyData();surface.DeepCopy(reader.GetOutput());v=vtk_to_numpy(surface.GetPoints().GetData()).astype(float);v=(v-origin)/dx;surface.GetPoints().SetData(numpy_to_vtk(v,deep=True));surface.Modified()
edges=vtk.vtkFeatureEdges();edges.SetInputData(surface);edges.BoundaryEdgesOn();edges.NonManifoldEdgesOn();edges.FeatureEdgesOff();edges.ManifoldEdgesOff();edges.Update();assert edges.GetOutput().GetNumberOfCells()==0,'Frozen closed STL is not closed/manifold'
sdf=vtk.vtkImplicitPolyDataDistance();sdf.SetInput(surface)
closed=np.fromfile(B/'frozen_inputs/step2/diagnostics/closed_flag_matrix.u8',dtype=np.uint8).reshape(nz,ny,nx)
zyx=np.argwhere(closed);centers=zyx[:,::-1].astype(float);coarse=centers[np.all((centers.astype(int)%2)==0,axis=1)]
# Signed metric cross-check on deterministic existing interior lattice points and a known exterior point.
interior=np.fromfile(B/'frozen_inputs/step2/diagnostics/palabos_native_flags.u8',dtype=np.uint8).reshape(nz,ny,nx)==3
ix=np.argwhere(interior)[::max(1,int(interior.sum())//128),::-1];signs=np.array([sdf.EvaluateFunction(p) for p in ix]);assert np.all(signs<0) and sdf.EvaluateFunction([-10.,-10.,-10.])>0,'SDF sign does not agree with frozen lumen'
def distances(p):return -np.array([sdf.EvaluateFunction(x) for x in p])
evals,evecs=np.linalg.eigh((verts.T@verts)/len(verts));normal=evecs[:,0]
def rot(axis,angle):
 axis=np.asarray(axis,dtype=float);axis/=np.linalg.norm(axis);x,y,z=axis;K=np.array([[0,-z,y],[z,0,-x],[-y,x,0]]);return np.eye(3)+math.sin(angle)*K+(1-math.cos(angle))*(K@K)
def align(a,b):
 cross=np.cross(a,b);n=np.linalg.norm(cross);d=np.dot(a,b)
 if n<1e-10:return np.eye(3) if d>0 else rot(np.eye(3)[np.argmin(np.abs(a))],math.pi)
 return rot(cross/n,math.atan2(n,d))
# Deterministic extremal vertices plus a uniform index sample for a cheap rejection pass.
dirs=np.array([[x,y,z] for x in (-1,0,1) for y in (-1,0,1) for z in (-1,0,1) if x or y or z]);subset=np.unique(np.r_[np.argmax(verts@dirs.T,axis=0),np.linspace(0,len(verts)-1,32,dtype=int)])
normals=[np.array([1.,0,0]),np.array([0.,1,0]),np.array([0.,0,1])]
for i in range(policy['normal_hemisphere_count']):
 z=(i+.5)/policy['normal_hemisphere_count'];angle=i*math.pi*(3-math.sqrt(5));normals.append(np.array([math.sqrt(1-z*z)*math.cos(angle),math.sqrt(1-z*z)*math.sin(angle),z]))
poses=[];broad_survivors=[];tested=0
for ni,n in enumerate(normals):
 for roll in policy['roll_degrees']:
  Q=rot(n,math.radians(roll))@align(normal,n);offset=verts[subset]@Q.T;points=np.rint(coarse[:,None,:]+offset[None,:,:]).astype(np.int32)
  valid=np.all((points>=0)&(points<np.array([nx,ny,nz])),axis=2);np.clip(points,0,np.array([nx-1,ny-1,nz-1]),out=points)
  inside=closed[points[:,:,2],points[:,:,1],points[:,:,0]].astype(bool)&valid;score=inside.sum(axis=1);tested+=len(coarse)
  top=np.argsort(score,kind='stable')[-3:]
  for j in top:poses.append((int(score[j]),int(ni),int(roll),coarse[j].copy(),Q.copy()))
  survivors=np.flatnonzero(score==len(subset))
  for j in survivors:
   points=np.rint(coarse[j]+verts@Q.T).astype(int)
   if np.all((points>=0)&(points<np.array([nx,ny,nz]))) and np.all(closed[points[:,2],points[:,1],points[:,0]]):broad_survivors.append((int(score[j]),int(ni),int(roll),coarse[j].copy(),Q.copy()))
 if ni%16==0:print('GEOMETRY_SEARCH',ni,len(normals),'poses',tested,'full_vertex_voxel_candidates',len(broad_survivors),flush=True)
ranked=sorted(poses,key=lambda x:(-x[0],x[1],x[2],tuple(x[3])));pool=broad_survivors+ranked[:policy['refinement_candidates']];seen=set();candidates=[]
for score,ni,roll,p,Q in pool:
 key=(tuple(p),ni,roll)
 if key in seen:continue
 seen.add(key);d=distances(p+verts@Q.T);candidates.append(dict(coarse_inside_count=score,normal_index=ni,roll_degrees=roll,center=p,rotation=Q,min_distance=float(d.min()),outside_vertices=int(np.sum(d<0))))
candidates.sort(key=lambda x:-x['min_distance'])
# Predeclared rigid-pose local refinement. It changes placement only before any solver run.
for ci,cand in enumerate(candidates[:policy['refinement_candidates']]):
 p=cand['center'].copy();Q=cand['rotation'].copy();best=cand['min_distance'];step=1.;angle=math.radians(5.)
 for iteration in range(policy['refinement_rounds']):
  trials=[]
  for ax in np.eye(3):
   for sign in (-1,1):trials.append((p+sign*step*ax,Q));trials.append((p,rot(ax,sign*angle)@Q))
  improved=False
  for pp,QQ in trials:
   val=float(distances(pp+verts@QQ.T).min())
   if val>best:p,Q,best=pp,QQ,val;improved=True
  if not improved:step*=.5;angle*=.5
 cand.update(center=p,rotation=Q,min_distance=best,outside_vertices=int(np.sum(distances(p+verts@Q.T)<0)))
 if ci%8==0:print('REFINE',ci,'best_min_clearance_LU',best,flush=True)
candidates.sort(key=lambda x:-x['min_distance'])
def polydata(points):
 poly=vtk.vtkPolyData();pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(points,deep=True));poly.SetPoints(pts);cells=vtk.vtkCellArray();packed=np.column_stack([np.full(len(faces),3),faces]).astype(np.int64);cells.SetCells(len(faces),numpy_to_vtkIdTypeArray(packed.ravel(),deep=True));poly.SetPolys(cells);return poly
def collision(poly):
 f=vtk.vtkCollisionDetectionFilter();f.SetInputData(0,poly);f.SetInputData(1,surface);t=vtk.vtkTransform();f.SetTransform(0,t);f.SetTransform(1,t);f.SetBoxTolerance(0);f.SetCellTolerance(0);f.SetCollisionModeToAllContacts();f.GenerateScalarsOff();f.Update();return f.GetNumberOfContacts()
verified=[]
for i,cand in enumerate(candidates[:64]):
 pts=cand['center']+verts@cand['rotation'].T;poly=polydata(pts);contacts=collision(poly);cand['wall_triangle_contacts']=int(contacts)
 if cand['min_distance']>=policy['CASE_A_C_min_wall_clearance_lu'] and contacts==0:verified.append(cand)
 if i==0:
  for name,obj in [('BEST_RIGID_PLACEMENT.vtp',poly),('FROZEN_LUMEN_LU.vtp',surface)]:
   w=vtk.vtkXMLPolyDataWriter();w.SetFileName(str(G/name));w.SetInputData(obj);w.Write()
 for key in ['center','rotation']:cand[key]=cand[key].tolist()
 # Prevent counting the same physical location multiple times in the small pack.
# Select a center case and 2-4 mutually separated safe placements; force center separation above full mesh bounding diameter.
pack=[];diam=2*float(np.linalg.norm(verts,axis=1).max())
for cand in verified:
 if all(np.linalg.norm(np.array(cand['center'])-np.array(other['center']))>diam+4 for other in pack):pack.append(cand)
 if len(pack)==4:break
best=candidates[0]
summary=dict(status='PASS' if verified and len(pack)>=2 else 'FAIL',RBC_GEOMETRIC_FIT='PASS' if verified else 'FAIL',CASE_A_fit='PASS' if verified else 'FAIL',CASE_C_fit='PASS' if len(pack)>=2 else 'FAIL',mesh=dict(vertices=len(verts),triangles=len(faces),nominal_scale=meta['scale'],nominal_diameter_um=2*meta['radius_m']*1e6,actual_volume_um3=vol,actual_area_um2=area,axis_extents_um=(np.ptp(verts,axis=0)*dxum).tolist(),maximum_vertex_extent_um=diam*dxum,principal_axis_extent_um=(np.ptp(verts@evecs,axis=0)*dxum).tolist()),policy=policy,SDF_sign_check='PASS',closed_surface_manifold='PASS',coarse_centers=len(coarse),coarse_poses=tested,coarse_vertex_count=len(subset),full_vertex_voxel_candidates=len(broad_survivors),continuous_verified_safe_placements=len(verified),best_pose={k:v for k,v in best.items()},CASE_A_spawn=verified[0] if verified else None,CASE_C_spawns=pack,small_pack_count=len(pack),search_seconds=time.monotonic()-t0,fluid_timesteps=0,RBC_timesteps=0,limitations=policy['search_scope'])
# Convert NumPy values without losing scientific evidence.
def convert(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,dict):return {k:convert(v) for k,v in x.items()}
 if isinstance(x,list):return [convert(v) for v in x]
 return x
summary=convert(summary);(G/'RBC_GEOMETRY_FIT_AUDIT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');(G/'PLACEMENT_CANDIDATES.json').write_text(json.dumps(convert(candidates[:64]),indent=2)+'\n');print(json.dumps(summary,indent=2))

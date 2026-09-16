"""First static geometry audit: C++ BVH against independent VTK and PCA."""
from pathlib import Path
import sys,json,ctypes,csv,time,hashlib
import numpy as np,h5py,vtk
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text());W=Path(P['local_work'])
sys.path.insert(0,str(R/'src'));from reference_wall_hydrodynamics_v0 import Geometry,frame
ptr=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');lib=ctypes.CDLL(str(W/'libwall_math.so'));lib.wall_geometry_open.argtypes=[ctypes.c_char_p];lib.wall_geometry_open.restype=ctypes.c_void_p;lib.wall_geometry_query.argtypes=[ctypes.c_void_p,ctypes.c_int,ptr,ptr,ptr]
def cpp(handle,x,a):
 out=np.empty((len(x),27));assert lib.wall_geometry_query(handle,len(x),np.ascontiguousarray(x),np.ascontiguousarray(a),out)==0;assert np.all(out[:,19]==1);return out
def save(name,obj):(R/name).write_text(json.dumps(obj,indent=2)+'\n')
assert not (R/'raw/REAL_GEOMETRY_FIRST.npz').exists(),'Preserve first static sample'
plan={'status':'FROZEN_BEFORE_STATIC_QUERY','seed':2026091602,'real_safe_points':10000,'curved_points':10000,
 'PCA_patch':'All triangles intersecting ball centered at closest point, radius2a; full vertex samples area/3; sorted original triangle IDs',
 'RMS_over_a_gate':.10,'P95_area_weighted_normal_spread_deg_gate':20.,'nearest_normal_angle_ambiguity_deg':20.,
 'triangle_normal_orientation':'raw geometric winding retained; flip only for fluid-facing diagnostic; gap direction n=(X-P)/distance',
 'curvature_estimator':'UNVERIFIED; do not infer curvature from RMS alone',
 'distance_gate_m':1e-12,'frame_identity_gate':1e-12,'PCA_numeric_rms_tolerance':1e-7,'PCA_normal_angle_tolerance_deg':1e-4}
save('contracts/GEOMETRY_AUDIT_PLAN.json',plan)
rng=np.random.default_rng(plan['seed']);a_values=np.array(list(json.loads((R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json').read_text())['radii_m'].values()))
geo=Geometry(R/'provenance/closed_geometry_m.stl');handle=lib.wall_geometry_open(str(R/'provenance/closed_geometry_m.stl').encode());assert handle
with h5py.File(R/'fields/FROZEN_FLOW_FIELD_V0.h5') as f:
 dims=f['dims'][:];origin=f['origin_m'][:];dx=float(f['dx_m'][()]);ids=f['linear_index'][:]
allij=np.column_stack([ids%dims[0],ids//dims[0]%dims[1],ids//(dims[0]*dims[1])]);centers=origin+dx*allij
order=rng.permutation(len(centers));positions=[];radii=[];states=[];sourceids=[];near_candidates=[]
for start in range(0,len(order),2000):
 ix=order[start:start+2000];x=centers[ix]+rng.uniform(-.2,.2,(len(ix),3))*dx;a=a_values[np.arange(start,start+len(ix))%3];out=cpp(handle,x,a)
 for j in np.where(out[:,13]>=0)[0]:
  if not geo.inside(x[j]):continue
  positions.append(x[j]);radii.append(a[j]);states.append(out[j]);sourceids.append(int(ids[ix[j]]))
  if len(positions)==10000:break
 print('real_safe',len(positions),'scanned',start+len(ix),flush=True)
 if len(positions)==10000:break
assert len(positions)==10000,'Insufficient safe interior centers; no radius changes'
x=np.array(positions);a=np.array(radii);out=np.array(states);np.savez_compressed(R/'raw/REAL_GEOMETRY_FIRST.npz',positions=x,radii=a,cpp_query=out,source_cell_ids=sourceids)
records=[];maxdist=0;maxpoint=0;maxrms=0;maxangle=0;class_disagreement=0
for i,(xx,aa,c) in enumerate(zip(x,a,out)):
 q=geo.query(xx,aa);de=abs(q['distance']-c[12]);pe=np.linalg.norm(q['closest']-c[:3]);re=abs(q['rms_over_a']-c[15]);ae=abs(q['normal_spread_deg']-c[16]);maxdist=max(maxdist,de);maxpoint=max(maxpoint,pe);maxrms=max(maxrms,re);maxangle=max(maxangle,ae);class_disagreement+=int(q['planar']!=bool(c[20]))
 row={'x_m':xx[0],'y_m':xx[1],'z_m':xx[2],'radius_m':aa,'triangle_id':int(c[14]),'gap_m':c[13],'epsilon':c[13]/aa,'rms_over_a':c[15],'normal_spread_deg':c[16],'normal_x':c[3],'normal_y':c[4],'normal_z':c[5],'nearest_normal_angle_deg':c[17],'patch_triangles':int(c[18]),'LOCAL_PLANE_VALID':'YES' if c[20] else 'NO','vtk_distance_error_m':de,'vtk_closest_point_error_m':pe,'independent_rms_error':re,'independent_p95_error_deg':ae,'source_cell_id':sourceids[i]};records.append(row)
 if (i+1)%1000==0:print('VTK_real',i+1,flush=True)
with (R/'validation/WALL_LOCAL_PLANE_VALIDITY.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
save('validation/REAL_GEOMETRY_AUDIT.json',{'status':'PASS' if maxdist<=1e-12 and maxrms<=1e-7 and maxangle<=1e-4 and class_disagreement==0 else 'FAIL','count':len(x),'safe_interior_verified_by_VTK':True,'max_distance_error_m':maxdist,'max_closest_point_error_m':maxpoint,'max_rms_error':maxrms,'max_normal_spread_error_deg':maxangle,'validity_class_disagreements':class_disagreement,'local_plane_valid_points':int(out[:,20].sum()),'local_plane_invalid_points':int((out[:,20]==0).sum()),'spatial_bounds_m':[x.min(axis=0).tolist(),x.max(axis=0).tolist()]})
# The original Case J population is never relocated or resized.
old=json.loads((R/'provenance/PREVIOUS_REAL_8_CASE_CONTRACT.json').read_text());jx=np.array(old['initial_positions_m']);ja=np.array(old['radii_m']);jo=cpp(handle,jx,ja)
jchecks=[]
for i in range(len(jx)):
 q=geo.query(jx[i],ja[i]);jchecks.append({'id':old['ids'][i],'cpp':jo[i].tolist(),'vtk_independent':{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in q.items()}})
save('validation/CASE_J_INITIAL_GEOMETRY.json',{'status':'PASS' if np.all(jo[:,20]) else 'BLOCKED_LOCAL_PLANE_VALIDITY','checks':jchecks,'population_changed':False})

# Select one genuine d50 near-wall point using only static criteria. Project each
# frozen flow node toward its nearest surface to epsilon0.1, then require an
# entirely valid existing flow interpolation cell and both geometry checks.
field_ids=set(map(int,ids));a50=a_values[1];candidates=[]
for start in range(0,len(centers),2000):
 pp=centers[start:start+2000];aa=np.full(len(pp),a50);oldq=cpp(handle,pp,aa);xx=oldq[:,:3]+oldq[:,3:6]*(1.1*a50);qq=cpp(handle,xx,aa)
 for j in np.where(qq[:,20]==1)[0]:
  ij=np.floor((xx[j]-origin)/dx).astype(int)
  if np.any(ij<0) or np.any(ij+1>=dims):continue
  corners=[int(ij[0]+u+dims[0]*(ij[1]+v+dims[1]*(ij[2]+w))) for u in [0,1] for v in [0,1] for w in [0,1]]
  if not all(k in field_ids for k in corners):continue
  independent=geo.query(xx[j],a50)
  if not independent['planar'] or not geo.inside(xx[j]):continue
  candidates.append({'position_m':xx[j].tolist(),'radius_m':a50,'epsilon':float(qq[j,13]/a50),'source_node_id':int(ids[start+j]),'rms_over_a':float(qq[j,15]),'normal_spread_deg':float(qq[j,16]),'triangle_id':int(qq[j,14])})
 if candidates:break
 print('single_start_candidates_scanned',start+len(pp),flush=True)
chosen=sorted(candidates,key=lambda c:(c['normal_spread_deg'],c['rms_over_a'],c['source_node_id']))[0] if candidates else None
save('contracts/CASE_I_INITIAL_SELECTION.json',{'status':'PASS' if chosen else 'BLOCKED_NO_VALID_NEAR_WALL_START','chosen':chosen,'selection':'First source-node batch with valid candidates; minimize normal spread, then RMS, then source node ID. Exact d50 radius and target epsilon0.1.', 'flow_stencil_all_8_corners_valid':bool(chosen),'number_candidates_in_selected_batch':len(candidates)})

radius=20e-6;sphere=vtk.vtkSphereSource();sphere.SetRadius(radius);sphere.SetThetaResolution(256);sphere.SetPhiResolution(128);sphere.Update();writer=vtk.vtkSTLWriter();writer.SetFileName(str(R/'provenance/SYNTHETIC_SPHERE_R20UM.stl'));writer.SetFileTypeToBinary();writer.SetInputData(sphere.GetOutput());writer.Write()
sg=Geometry(R/'provenance/SYNTHETIC_SPHERE_R20UM.stl');sh=lib.wall_geometry_open(str(R/'provenance/SYNTHETIC_SPHERE_R20UM.stl').encode());assert sh
normals=rng.normal(size=(10000,3));normals/=np.linalg.norm(normals,axis=1)[:,None];sa=np.resize(a_values,10000);h=sa*np.exp(rng.uniform(np.log(.01),np.log(5),10000));sx=normals*(radius-sa-h)[:,None];so=cpp(sh,sx,sa);analytic_distance=radius-np.linalg.norm(sx,axis=1);distance_error=[];point_error=[];analytic_error=abs(so[:,12]-analytic_distance)
for i in range(10000):
 p,cell,dist=sg.closest(sx[i]);distance_error.append(abs(dist-so[i,12]));point_error.append(np.linalg.norm(p-so[i,:3]))
np.savez_compressed(R/'raw/CURVED_GEOMETRY_FIRST.npz',positions=sx,radii=sa,cpp_query=so,analytic_distance=analytic_distance)
bound=radius*(1-np.cos(np.pi/127))*2+5e-12
save('validation/CURVED_GEOMETRY_AUDIT.json',{'status':'PASS' if max(distance_error)<=1e-12 and max(analytic_error)<=bound else 'FAIL','count':10000,'analytic_sphere_radius_m':radius,'max_CPP_VTK_distance_error_m':max(distance_error),'max_CPP_VTK_closest_point_error_m':max(point_error),'max_STL_analytic_distance_error_m':float(max(analytic_error)),'mesh_discretization_distance_bound_m':bound,'analytic_distance_gate':'Mesh sagitta bound, distinct from 1e-12 CPP/VTK same-mesh identity','all_local_frames_orthonormal':'checked in rotation audit','local_plane_valid_points':int(so[:,20].sum())})
print('STATIC_GEOMETRY_COMPLETE',flush=True)

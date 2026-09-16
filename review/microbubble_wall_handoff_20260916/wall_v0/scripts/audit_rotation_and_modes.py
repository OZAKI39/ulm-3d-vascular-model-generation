from pathlib import Path
import ctypes,json,sys,csv
import numpy as np,h5py
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text());W=Path(P['local_work']);sys.path.insert(0,str(R/'src'))
from reference_wall_hydrodynamics_v0 import Lookup,frame
ptr=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');lib=ctypes.CDLL(str(W/'libwall_math.so'));lib.wall_lookup_open.argtypes=[ctypes.c_char_p];lib.wall_lookup_open.restype=ctypes.c_void_p
table=R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5';handle=lib.wall_lookup_open(str(table).encode());assert handle
lib.wall_rotation_audit.argtypes=[ctypes.c_int,ptr,ptr,ptr,ptr,ptr,ptr];lib.wall_single_response.argtypes=[ctypes.c_void_p,ctypes.c_int,ptr,ptr]
rng=np.random.default_rng(2026091603);n=10000;normals=rng.normal(size=(n,3));normals/=np.linalg.norm(normals,axis=1)[:,None];ee=np.exp(rng.uniform(np.log(.001),np.log(5),n));M=np.ascontiguousarray(Lookup(table).excess(ee)+np.eye(6));v=rng.normal(size=(n,6));frames=np.zeros((n,3,3));f1=np.zeros((n,6));f2=f1.copy();assert lib.wall_rotation_audit(n,normals,M,v,frames,f1,f2)==0
rel=np.linalg.norm(f1-f2,axis=1)/np.linalg.norm(f2,axis=1);ortho=np.max(abs(frames.transpose(0,2,1)@frames-np.eye(3)));deterr=np.max(abs(np.linalg.det(frames)-1));pyframes=np.array([frame(x) for x in normals]);frameerr=np.max(abs(frames-pyframes));assert np.max(rel)<=1e-12 and ortho<=1e-12 and deterr<=1e-12 and frameerr<=1e-12
np.savez_compressed(R/'raw/ROTATION_FIRST.npz',normals=normals,matrices=M,velocities=v,frames=frames,force_matrix=f1,force_two_rotations=f2)
status={'status':'PASS','normal_count':n,'max_relative_rotation_identity_error':float(max(rel)),'max_orthonormal_error':float(ortho),'max_determinant_error':float(deterr),'max_independent_frame_error':float(frameerr),'no_reflection':True}
(R/'validation/ROTATION_AUDIT.json').write_text(json.dumps(status,indent=2)+'\n')
inputs=[];refs=[];info=[]
with h5py.File(R/'provenance/reference_audit/WALL_REFERENCE_MOBILITY_MATRICES.h5') as f:
 for i in range(len(f['epsilon'])):
  if f['implementation'].asstr()[i]!='RMBW_LUBRICATION' or f['viscosity_Pa_s'][i]!=.001:continue
  e=float(f['epsilon'][i]);a=float(f['radius_m'][i]);size=f['size'].asstr()[i]
  if e not in [.001,.005,.01,.05,.1,.2,.5,1.,5.]:continue
  for case,axis in [('B_normal',2),('C_tangential',0),('D_TR',4),('E_parallel_rotation',3),('E_normal_rotation',5)]:
   load=np.zeros(6);load[axis]=(6*np.pi*.001*a*1e-4 if axis<3 else 8*np.pi*.001*a**3*10)
   inputs.append([a,e*a,*load,0,0,1]);refs.append(f['matrix_SI'][i]@load);info.append((case,size,e,a))
inp=np.array(inputs);ref=np.array(refs);vel=np.empty_like(ref);assert lib.wall_single_response(handle,len(inp),inp,vel)==0
records=[]
for i,(case,size,e,a) in enumerate(info):
 scale=np.array([1.]*3+[a]*3);err=np.linalg.norm((vel[i]-ref[i])*scale)/np.linalg.norm(ref[i]*scale);components=np.where(abs(ref[i])>1e-30,abs(vel[i]-ref[i])/np.maximum(abs(ref[i]),1e-300),abs(vel[i]-ref[i]));records.append({'case':case,'size':size,'epsilon':e,'scaled_velocity_error':err,'max_nonzero_component_error':float(max(components)),'status':'PASS' if max(err,max(components))<=1e-3 else 'FAIL',**{'computed_'+str(j):vel[i,j] for j in range(6)},**{'reference_'+str(j):ref[i,j] for j in range(6)}})
with (R/'validation/FLAT_WALL_MODE_RESPONSES.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
np.savez_compressed(R/'raw/FLAT_MODE_SOLVER_FIRST.npz',inputs=inp,cpp_velocity=vel,RMBW_reference_velocity=ref)
summary={'status':'PASS' if all(x['status']=='PASS' for x in records) else 'FAIL','cases':len(records),'maximum_scaled_velocity_error':max(x['scaled_velocity_error'] for x in records),'maximum_component_error':max(x['max_nonzero_component_error'] for x in records),'solver':'actual production rigid::solve with isolated bulk plus wall excess','BULK_DOUBLE_COUNTING':'NO','rotation_qualification':'PASS_WITH_REFERENCE_LIMITATION'}
(R/'validation/FLAT_WALL_MODE_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));assert summary['status']=='PASS'

"""Additional static search only; preserve first failed source-node search."""
from pathlib import Path
import json,sys,ctypes
import numpy as np,h5py
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text());W=Path(P['local_work']);sys.path.insert(0,str(R/'src'))
from reference_wall_hydrodynamics_v0 import Geometry
g=Geometry(R/'provenance/closed_geometry_m.stl');a=json.loads((R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json').read_text())['radii_m']['d50']
plan={'status':'FROZEN_BEFORE_SUPPLEMENTAL_STATIC_SEARCH','sphere_radius_m':a,'surface_points':'all original triangle centroids','epsilon_candidates':[.05,.2],'normal_sign_candidates':[-1,1],'geometry_gates_unchanged':True,'initial_failed_search_preserved':True}
(R/'contracts/CASE_I_SUPPLEMENTAL_SEARCH_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
lib=ctypes.CDLL(str(W/'libwall_math.so'));lib.wall_geometry_open.argtypes=[ctypes.c_char_p];lib.wall_geometry_open.restype=ctypes.c_void_p;h=lib.wall_geometry_open(str(R/'provenance/closed_geometry_m.stl').encode());ptr=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');lib.wall_geometry_query.argtypes=[ctypes.c_void_p,ctypes.c_int,ptr,ptr,ptr]
with h5py.File(R/'fields/FROZEN_FLOW_FIELD_V0.h5') as f:dims=f['dims'][:];origin=f['origin_m'][:];dx=float(f['dx_m'][()]);ids=set(map(int,f['linear_index'][:]))
found=[];tested=0;geomvalid=0
for e in [.05,.2]:
 for sign in [-1,1]:
  for start in range(0,len(g.tri),2000):
   xx=np.ascontiguousarray(g.centers[start:start+2000]+sign*a*(1+e)*g.normal[start:start+2000]);aa=np.full(len(xx),a);out=np.empty((len(xx),27));assert lib.wall_geometry_query(h,len(xx),xx,aa,out)==0;tested+=len(xx)
   for j in np.where((out[:,19]==1)&(out[:,20]==1)&(out[:,13]/a>=.05-1e-10)&(out[:,13]/a<=.2+1e-10))[0]:
    geomvalid+=1;ij=np.floor((xx[j]-origin)/dx).astype(int)
    if np.any(ij<0) or np.any(ij+1>=dims):continue
    corners=[int(ij[0]+u+dims[0]*(ij[1]+v+dims[1]*(ij[2]+w))) for u in [0,1] for v in [0,1] for w in [0,1]]
    if not all(k in ids for k in corners) or not g.inside(xx[j]):continue
    q=g.query(xx[j],a)
    if not q['planar']:continue
    found.append({'position_m':xx[j].tolist(),'radius_m':a,'epsilon':out[j,13]/a,'triangle_id':int(out[j,14]),'rms_over_a':out[j,15],'normal_spread_deg':out[j,16]})
  print('searched',e,sign,'total',tested,'geometry_valid',geomvalid,'eligible',len(found),flush=True)
chosen=sorted(found,key=lambda x:(x['normal_spread_deg'],x['rms_over_a'],x['triangle_id']))[0] if found else None
(R/'contracts/CASE_I_SUPPLEMENTAL_SELECTION.json').write_text(json.dumps({'status':'PASS' if chosen else 'BLOCKED_NO_VALID_NEAR_WALL_START','chosen':chosen,'tested':tested,'geometry_valid_candidates':geomvalid,'eligible_candidates':len(found),'not_an_exhaustive_continuous_domain_proof':True},indent=2)+'\n')

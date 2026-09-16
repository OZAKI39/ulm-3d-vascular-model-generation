from pathlib import Path
import sys,json,ctypes,numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'));from reference_wall_hydrodynamics_v0 import Geometry
W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work'])
seed=2026091604;plan={'seed':seed,'N':200,'radius_m':.9683592065545495e-6,'synthetic_sphere_radius_m':20e-6,'families':['interior safe segments','outside safe same-side segments','endpoints-safe through-sphere crossings','endpoint overlap'],'purpose':'continuous swept sphere/triangle adversarial audit; no simulation timestep'};(R/'contracts/SWEPT_STATIC_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
rng=np.random.default_rng(seed);a=np.full(200,plan['radius_m']);directions=rng.normal(size=(200,3));directions/=np.linalg.norm(directions,axis=1)[:,None];x=np.zeros((200,3));y=x.copy()
x[:50]=directions[:50]*10e-6;y[:50]=directions[50:100]*10e-6
x[50:100]=directions[50:100]*25e-6;y[50:100]=directions[50:100]*28e-6
x[100:150]=directions[100:150]*25e-6;y[100:150]=-directions[100:150]*25e-6
x[150:]=directions[150:]*20e-6;y[150:]=directions[150:]*25e-6
lib=ctypes.CDLL(str(W/'libwall_math.so'));lib.wall_geometry_open.argtypes=[ctypes.c_char_p];lib.wall_geometry_open.restype=ctypes.c_void_p;handle=lib.wall_geometry_open(str(R/'provenance/SYNTHETIC_SPHERE_R20UM.stl').encode());assert handle
ptr=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');ip=np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS');lib.wall_geometry_swept.argtypes=[ctypes.c_void_p,ctypes.c_int,ptr,ptr,ptr,ip];out=np.zeros(200,np.int32);assert lib.wall_geometry_swept(handle,200,x,y,a,out)==0
geo=Geometry(R/'provenance/SYNTHETIC_SPHERE_R20UM.stl');ref=np.array([geo.swept_safe(xx,yy,aa) for xx,yy,aa in zip(x,y,a)]);expected=np.r_[np.ones(100,bool),np.zeros(100,bool)];assert np.array_equal(ref,out) and np.array_equal(ref,expected)
np.savez_compressed(R/'raw/SWEPT_STATIC_FIRST.npz',x=x,y=y,radii=a,cpp_safe=out,VTK_safe=ref,expected=expected)
(R/'validation/SWEPT_STATIC_AUDIT.json').write_text(json.dumps({'status':'PASS','N':200,'safe':100,'rejected':100,'endpoints_safe_but_interior_crossing_rejections':50,'independent_VTK_comparison':'PASS','same_frozen_production_segment_algorithm':True},indent=2)+'\n');print('SWEPT_STATIC_PASS')

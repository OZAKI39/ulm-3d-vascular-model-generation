from pathlib import Path
import sys,json,time,ctypes
import numpy as np
sys.path.insert(0,str(Path.cwd()/'particle_3d/src'))
from particle_3d.particle82_tracers import init_environment
from particle_3d.particle82_point_native import NativePointTracer
E=init_environment();solver=NativePointTracer(E);base=Path.cwd().parent
p=np.load(base/'native_point_smoke/birth_positions.npz')['positions'][21]
checks=[]
for eps in [1e-11,1e-12,1e-13]:
 row=solver.trace(p,error=eps,step_m=2e-7);a=row.pop('path');last=E.field.sample(a[-1,1:]);_,distance=E.wall.nearest_center_triangle(a[-1,1:])
 row.update(error=eps,elapsed_s=float(a[-1,0]),saved_states=len(a),position_m=a[-1,1:].tolist(),velocity_m_s=last.velocity_m_s.tolist(),tetra_id=last.tetra_id,velocity_gradient_s_inv=last.velocity_gradient_s_inv.tolist(),divergence_s_inv=float(np.trace(last.velocity_gradient_s_inv)),distance_to_wall_m=distance,nearby_last_motion_m=float(np.linalg.norm(a[-1,1:]-a[-100,1:])))
 checks.append(row);print(json.dumps(row),flush=True)
# Independent field sampler parity at 1024 original tetra barycenters.
ids=np.linspace(0,len(E.field.tetra)-1,1024,dtype=int);errors=[]
for cell in ids:
 x=np.ascontiguousarray(E.field.points_m[E.field.tetra[cell]].mean(axis=0));v=np.zeros(3)
 nativecell=solver.lib.p82_sample(x.ctypes.data,int(cell),v.ctypes.data);s=E.field.sample(x)
 assert nativecell==s.tetra_id
 errors.append(float(np.max(np.abs(v-s.velocity_m_s))))
record=dict(stagnant_seed_21=checks,field_sample_count=len(ids),velocity_max_abs_difference_m_s=max(errors),field_absolute_tolerance_m_s=1e-16,field_sampler_parity_pass=max(errors)<1e-16)
(base/'POINT_TRACER_IMPLEMENTATION_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n')
print('FIELD_PARITY',max(errors),flush=True)

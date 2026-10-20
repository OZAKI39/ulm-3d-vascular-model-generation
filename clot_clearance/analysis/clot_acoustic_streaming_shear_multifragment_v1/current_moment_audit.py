"""Actual moment at saved accepted states, distinct from nominal reference torque."""
from common import *
from shear_loading import ShearLoad
import csv
ref=load_npz(O/'provenance/reference.npz');records=[];errors=[]
for name in ['stokes','inverse']:
 load=ShearLoad(ref,name);seen=set()
 for path in sorted((O/'cases'/name/'keys').glob('*.npz')):
  if 'failed' in path.name:continue
  s=load_npz(path);step=int(s['step'])
  if step in seen:continue
  seen.add(step);force=s['external_force'];actual=np.cross(s['x'],force).sum(0);reference=np.cross(ref['X'],force).sum(0)
  virtual_face_centers=s['x'][load.tri].mean(1);virtual_face_moment=np.cross(virtual_face_centers,load.faces(float(s['time']))).sum(0)
  error=float(np.linalg.norm(actual-virtual_face_moment));assert error<1e-20;errors.append(error)
  row=dict(case=name,step=step,time_s=float(s['time']),adjoint_current_moment_error_Nm=error)
  for j,axis in enumerate('xyz'):row['actual_moment_'+axis+'_Nm']=actual[j];row['reference_moment_'+axis+'_Nm']=reference[j]
  records.append(row)
records.sort(key=lambda r:(r['case'],r['step']))
with (O/'CURRENT_MOMENT_KEY_STATES.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
dump(O/'MOMENT_CONVENTION_AUDIT.json',dict(status='PASS',key_states=len(records),max_current_moment_adjoint_error_Nm=max(errors),history_columns={'net_moment_z_Nm':'Nominal reference-coordinate torque sum X cross F; not the current-position torque','moment_impulse_z_Nms':'Time integral of the nominal reference-coordinate torque','impulse_x_Ns':'Actual prescribed-force impulse; position independent','absolute_impulse_Ns':'Integral of sum of nodal force magnitudes; does not hide counter-shear cancellation'},current_moment='Computed directly from actual x(t) cross native external force for every unique saved accepted FP64 state.',surface_caveat='Moved original triangle barycenters establish an algebraic moment/power identity only. Disconnected/deformed original triangles are not claimed to be actual newly exposed fluid boundaries.',current_angular_impulse='UNAVAILABLE_AT_EVERY_ACCEPTED_STEP; not reconstructed from insufficient sampled positions'))
print('CURRENT_MOMENT_AUDIT_PASS',len(records))

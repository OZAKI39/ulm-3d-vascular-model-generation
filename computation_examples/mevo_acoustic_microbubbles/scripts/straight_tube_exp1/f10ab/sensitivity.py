"""Central differences of invariant shape outputs, with step and solver noise bounds."""
import copy,h5py,numpy as np
from .common import *
from .metrics import measure
from .plan import task,freeze

def plan():
 admissions=read(BASE/'analysis/ADMISSION.json');valid={r['state_id'] for r in admissions if r['action']=='A7_X210' and r['NUMERICALLY_RESOLVED']}
 tasks=[];actions=read(BASE/'ACTIONS.json')
 for state in ['S2_r0','S4_r0']:
  if state not in valid:continue
  for parameter in ['amp','phase']:
   for h in [1,.5]:
    for sign in [-1,1]:
     rec=copy.deepcopy(actions['A7_X210']);seg=rec['segments'][0]
     if parameter=='amp':
      for i in [2,3]:seg['amplitude_Pa_equivalent'][i]*=1+sign*h*.05
     else:seg['phase_rad'][3]+=np.deg2rad(sign*h*6)
     name=f'J_{parameter}_{"p" if sign>0 else "m"}_{"full" if h==1 else "half"}'
     tasks.append(task('sensitivity',state,name,.002,schedule=rec))
 freeze('sensitivity',tasks,500,'Central ±5/2.5% common X-pair amplitude and ±6/3 degree relative phase around A7; S2/S4 discovery shapes; evaluate 0.5ms and 2ms from actual states. No global four-channel controllability claim.')

OUTPUTS=['COM X / 20nm','COM Y / 20nm','log(lambda1/lambda3) / .01','log(lambda2/lambda3) / .01','Rg fraction / .01','similarity RMS/Rg0 / .01']
def output(row,t):
 with h5py.File(BASE/row['trajectory_path']) as h:
  ts=h['time_s'][:];i=int(np.searchsorted(ts,t-1e-13));assert i<len(ts) and abs(ts[i]-t)<1e-13
  x0=h['positions_m'][0];x=h['positions_m'][i];a=measure(x0,x0);m=measure(x0,x)
 la=a['eigenvalues_3d_m2'];l=m['eigenvalues_3d_m2'];c=(m['center_3d_m']-a['center_3d_m'])/20e-9
 return np.r_[c[:2],(np.log(l[:2]/l[2])-np.log(la[:2]/la[2]))/.01,(m['Rg_3d_m']/a['Rg_3d_m']-1)/.01,m['similarity_rmsd_3d_m']/a['Rg_3d_m']/.01]

def analyze():
 rows=read(BASE/'sensitivity/INDEX.json');refs=read(BASE/'discovery/INDEX.json');refines=read(BASE/'refinement/INDEX.json');cases=[]
 for state in ['S2_r0','S4_r0']:
  group=[r for r in rows if r.get('state_id')==state]
  if len(group)!=8 or any(r['status']!='PASS' for r in group):
   cases.append(dict(state_id=state,status='INCONCLUSIVE',reason='Incomplete central-difference stencil'));continue
  base=next(r for r in refs if r['state_id']==state and r['action']=='A7_X210')
  finer=[r for r in refines if r.get('state_id')==state and r.get('action')=='A7_X210' and r['status']=='PASS']
  for t in [.0005,.002]:
   noise=np.max([abs(output(r,t)-output(base,t)) for r in finer],axis=0)
   controls=read(BASE/'sensitivity_checks/INDEX.json') if (BASE/'sensitivity_checks/COMPLETE.json').exists() else []
   matched=[r for r in controls if r.get('state_id')==state]
   stencil_checks_complete=len(matched)==8 and all(r['status']=='PASS' for r in matched)
   if stencil_checks_complete:
    for r in matched:
     ref=next(q for q in group if q['action']==r['action'])
     noise=np.maximum(noise,abs(output(r,t)-output(ref,t)))
   matrices=[]
   for suffix,h in [('full',1),('half',.5)]:
    cols=[]
    for parameter in ['amp','phase']:
     plus=next(r for r in group if r['action']==f'J_{parameter}_p_{suffix}');minus=next(r for r in group if r['action']==f'J_{parameter}_m_{suffix}')
     cols.append((output(plus,t)-output(minus,t))/(2*h))
    matrices.append(np.array(cols).T)
   j,k=matrices;difference=np.linalg.norm(j-k,2);noise_norm=5*np.linalg.norm(noise)*2
   singular=np.linalg.svd(k,compute_uv=False);shape=k[[2,3,5]];ss=np.linalg.svd(shape,compute_uv=False)
   threshold=max(difference,noise_norm,1e-6);finite_difference_converged=difference<=.2*max(np.linalg.norm(k,2),1e-12)
   # A Taylor direction must survive BOTH stencil and solver uncertainty.
   rank=int(np.count_nonzero(singular>threshold));srank=int(np.count_nonzero(ss>threshold))
   u,s,vt=np.linalg.svd(k[:2],full_matrices=True);candidate=vt[-1];shape_gain=float(np.linalg.norm(shape@candidate));drift=float(np.linalg.norm(k[:2]@candidate))
   cases.append(dict(state_id=state,time_s=t,status='PASS' if finite_difference_converged and stencil_checks_complete else 'INCONCLUSIVE',outputs=OUTPUTS,
    stencil_numerical_controls_complete=stencil_checks_complete,
    normalized_jacobian=k,coarse_difference_jacobian=j,normalization='input 1 = 5% pair amplitude or 6deg relative phase; outputs explicit per labels',
    numerical_noise=noise,stencil_discrepancy=float(difference),rank_threshold=float(threshold),effective_rank=rank,shape_effective_rank=srank,
    singular_values=singular,shape_singular_values=ss,condition_number=float(singular[0]/singular[-1]) if singular[-1]>threshold else None,
    minimum_translation_direction=candidate,shape_gain_along_direction=shape_gain,translation_gain_along_direction=drift,
    scope='Two-input local diagnostic at a fixed initial shape and time; no global nonlinear controllability proof; direction is not a validated control trajectory.'))
 write(BASE/'analysis/SENSITIVITY.json',dict(cases=cases,identity=identity(),training_ready=False))
 print([(r['state_id'],r.get('time_s'),r['status'],r.get('shape_effective_rank')) for r in cases])

def checks_plan():
 rows=read(BASE/'sensitivity/INDEX.json');tasks=[]
 for r in rows:
  if r['status']=='PASS' and r['action'].endswith('_half'):
   for variant in ['DT05','FINE']:
    tasks.append(task('sensitivity_checks',r['state_id'],r['action'],.002,variant,r['task']['schedule']))
 freeze('sensitivity_checks',tasks,450,'Both half-stencil endpoints receive original dt/cache controls; no rank conclusion from baseline noise alone.')

def verify_checks():
 from .analyze import numerical
 refs=read(BASE/'sensitivity/INDEX.json');checks=[]
 for b in read(BASE/'sensitivity_checks/INDEX.json'):
  a=next(r for r in refs if r['state_id']==b['state_id'] and r['action']==b['action'])
  c=numerical(a,b);c.update(state_id=a['state_id'],action=a['action'],variant=b['variant']);checks.append(c)
 write(BASE/'analysis/SENSITIVITY_NUMERICAL.json',checks)
 analyze()
 data=read(BASE/'analysis/SENSITIVITY.json')
 for case in data['cases']:
  numeric=all(c['original']['admitted'] for c in checks if c['state_id']==case['state_id'])
  case['original_full_window_gates_pass']=numeric
  if not numeric:case['status']='INCONCLUSIVE'
 write(BASE/'analysis/SENSITIVITY.json',data)
if __name__=='__main__':
 import sys
 plan() if len(sys.argv)>1 and sys.argv[1]=='plan' else analyze()

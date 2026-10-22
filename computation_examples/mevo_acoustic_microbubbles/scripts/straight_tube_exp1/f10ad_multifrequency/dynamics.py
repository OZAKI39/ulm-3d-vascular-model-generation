"""Limited multi-frequency force, complete trajectory and convergence admission."""
import time,h5py,numpy as np
from scipy.linalg import cho_factor,cho_solve
from mevo_mb.fast_sim.nearfield import pair_reference
from mevo_mb.fast_sim.accelerated_contact import resistance
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from scripts.straight_tube_exp1.f10ac_qd.common import *
from scripts.straight_tube_exp1.f10ac_qd.run import task,batch
from .interface import context,make_schedule,FREQUENCIES

def force_audit():
 rows=[]
 with h5py.File(BASE/'INITIAL_STATES.h5') as h:x0=h['S2_r0'][:]
 for f in FREQUENCIES:
  for shell in [1.,.8,1.2]:
   c=context(f,shell);basis=c['basis'];response=c['response'];phys=c['physics']
   for offset in [-5e-6,0,5e-6]:
    x=x0+[0,offset,0]
    for name,amps,phases in [('X180',[0,0,700,700],[0,0,0,np.pi]),('Y90',[700,700,0,0],[0,np.pi/2,0,0]),('FOUR',[700]*4,[0,np.pi/2,0,np.pi])]:
     try:
      e=AuditEngine(x,phys,response,basis,c['cavity'],policy=policy());e.control(make_schedule(amps,phases,.0005,f).segments[0].event());gpu=e.fields(x)
      pp,gg=basis.query_basis(x);u=basis.vector(amps,phases);p=u@pp;grad=np.einsum('i,inj->nj',u,gg);v=response.query(p)
      primary=-.5*np.real(v[:,0,None]*np.conj(grad));chi=response.query(np.array([complex(response.knots[0])]))[0,0]/response.knots[0];secondary=pair_reference(x,v,p,phys,chi)
      err1=float(abs(primary-gpu['primary_force_N']).max());err2=float(abs(secondary-gpu['secondary_force_N']).max());vel=cho_solve(cho_factor(resistance(x,R)['A']),(primary+secondary).ravel()/(6*np.pi*.001*R)).reshape(x.shape)
      e.gate_records=[];step=e.step(1e-6)
      r=dict(frequency=f,shell=shell,offset_y_m=offset,action=name,pressure_range_Pa=[float(abs(p).min()),float(abs(p).max())],primary_mean_N=primary.mean(0),primary_max_N=float(np.linalg.norm(primary,axis=1).max()),secondary_max_N=float(np.linalg.norm(secondary,axis=1).max()),free_velocity_mean_m_s=vel.mean(0),relative_velocity_rms_m_s=float(np.sqrt(np.mean((vel-vel.mean(0))**2))),constrained_velocity_mean_m_s=step['result']['velocity_m_s'].mean(0),single_step_gates=step['gates'],primary_cpu_gpu_error_N=err1,secondary_cpu_gpu_error_N=err2,PASS=err1<1e-18 and err2<1e-18)
      dst=AD/'static'/f'{int(f)}_{shell:.1f}_{offset:.2g}_{name}.npz';dst.parent.mkdir(exist_ok=True);np.savez_compressed(dst,positions_m=x,pressure_Pa=p,gradient=grad,harmonics_m3=v,primary_force_N=primary,secondary_force_N=secondary,free_velocity_m_s=vel,constrained_velocity_m_s=step['result']['velocity_m_s'])
     except Exception as ex:r=dict(frequency=f,shell=shell,offset_y_m=offset,action=name,PASS=False,reason=repr(ex))
     rows.append(r)
 write(AD/'FORCE_AUDIT.json',rows);return rows

def main():
 verify();start=time.perf_counter()
 from scripts.straight_tube_exp1.f10ac_qd.profile_fields import main as profile_fields
 profile_fields();forces=force_audit();plans=[]
 actions=[('OFF',[0]*4,[0]*4),('SINGLE_X',[0,0,700,0],[0]*4),('X180',[0,0,700,700],[0,0,0,np.pi]),('Y90',[700,700,0,0],[0,np.pi/2,0,0]),('FOUR',[700]*4,[0,np.pi/2,0,np.pi])]
 for f in [2400000.,2600000.]:
  if not read(AD/str(int(f))/'PHYSICS_GATE.json')['PASS']:continue
  for state in ['S0_r0','S2_r0']:
   for name,a,p in actions:
    plans.append(task(f'freq_{int(f)}_{state}_{name}',state,a,p,.0005,frequency=f))
 write(AD/'DYNAMICS_PREREGISTRATION.json',dict(short_plans=plans,extension='2ms X180/Y90/FOUR after complete short DT1/DT05/FINE original convergence, on S2; at all three frequencies',identity=identity()))
 rows=batch(plans,4);write(AD/'SHORT.json',rows);refs=[]
 for r in rows:
  if r['status']=='PASS' and (r['task']['state']=='S2_r0' and r['id'].split('_')[-1] in ['X180','Y90','FOUR'] or r['task']['state']=='S0_r0' and r['id'].endswith('X180')):
   for variant in ['DT05','FINE']:
    t=r['task'].copy();t.update(id=t['id']+'_'+variant,variant=variant);refs.append(t)
 refined=batch(refs,4);write(AD/'REFINED.json',refined);checks=[]
 for r in rows:
  cc=[]
  for b in refined:
   if b['id'].startswith(r['id']+'_'):
    result=compare((BASE/r['trajectory_path']).parent,(BASE/b['trajectory_path']).parent);result.update(reference=r['id'],comparison=b['id']);cc.append(result);checks.append(result)
  r['numerical_admitted']=len(cc)==2 and all(c['admitted'] for c in cc)
 write(AD/'SHORT_ADMISSION.json',rows);write(AD/'NUMERICAL_COMPARISONS.json',checks)
 long=[]
 for f in FREQUENCIES:
  for name,a,p in actions[2:]:
   admitted=f==2500000. or any(r['task'].get('frequency')==f and r['task']['state']=='S2_r0' and r['id'].endswith(name) and r['numerical_admitted'] for r in rows)
   if admitted:long.append(task(f'frequency_long_{int(f)}_{name}',amps=a,phase=p,frequency=f))
 results=batch(long,4);write(AD/'LONG.json',results)
 # Shell and Y near-node spatial robustness: actual dynamics, not force-only inference.
 stress=[]
 for r in results:
  if r['status']!='PASS':continue
  if r['id'].endswith('X180'):
   for shell in [.8,1.2]:t=r['task'].copy();t.update(id=t['id']+f'_shell{shell:.1f}',shell=shell);stress.append(t)
  if r['id'].endswith('Y90'):
   for offset in [-5e-6,5e-6]:t=r['task'].copy();t.update(id=t['id']+f'_offset{offset:.2g}',offset_m=[0,offset,0]);stress.append(t)
 # Every baseline long path has half-step and fine-cache checks.
 for r in results:
  if r['status']=='PASS':
   for v in ['DT05','FINE']:t=r['task'].copy();t.update(id=t['id']+'_'+v,variant=v);stress.append(t)
 rs=batch(stress,4);write(AD/'ROBUSTNESS.json',rs)
 comparisons=[]
 for r in results:
  for b in rs:
   if b['id'] in [r['id']+'_DT05',r['id']+'_FINE']:
    z=compare((BASE/r['trajectory_path']).parent,(BASE/b['trajectory_path']).parent);z.update(reference=r['id'],comparison=b['id']);comparisons.append(z)
 write(AD/'LONG_CONVERGENCE.json',comparisons)
 write(AD/'DYNAMICS_COMPLETE.json',dict(wall_s=time.perf_counter()-start,short=len(rows),refinements=len(refined),long=len(results),robustness=len(rs),identity=identity()))
if __name__=='__main__':main()

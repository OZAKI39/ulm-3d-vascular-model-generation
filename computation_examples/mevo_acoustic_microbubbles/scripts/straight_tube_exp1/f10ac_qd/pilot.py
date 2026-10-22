"""Independent environment regression and measured 1/2/4 worker pilot."""
import time,h5py,numpy as np,os,importlib.metadata as md
from .common import *
from .run import task,batch
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context,make_schedule
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from mevo_mb.control.scheduled_pzt import Schedule

def main():
 verify();start=time.perf_counter();records=[];parallel=[]
 # Exact historical trajectories rerun with same source amplitudes, time, cache and initial IDs.
 actions=read(AB/'ACTIONS.json')
 plan=[]
 for state,action in [('S0_r0','A7_X210'),('S2_r0','A5_X180')]:
  seg=actions[action]['segments'][0];plan.append(task('regression_'+state,state,seg['amplitude_Pa_equivalent'],seg['phase_rad'],.002))
 write(BASE/'plans/pilot.json',plan)
 rows=batch(plan,2);reg=[]
 historical=list(AB.glob('discovery/*/SUMMARY.json'))
 for r,(_,action) in zip(rows,[('S0_r0','A7_X210'),('S2_r0','A5_X180')]):
  ref=next(read(p) for p in historical if read(p)['state_id']==r['task']['state'] and read(p)['action']==action)
  with h5py.File(BASE/r['trajectory_path']) as a,h5py.File(AB/ref['trajectory_path']) as b:
   errors={k:float(np.max(abs(a[k][:]-b[k][:]),initial=0)) for k in ['time_s','positions_m','force_samples/pressure_Pa','force_samples/primary_force_N','force_samples/secondary_force_N']}
   same_contacts=bool(np.array_equal(a['force_samples/pair_contact_reaction_N'][:],b['force_samples/pair_contact_reaction_N'][:]))
  reg.append(dict(id=r['id'],errors=errors,contact_reactions_bitwise_equal=same_contacts,PASS=r['status']=='PASS' and errors['positions_m']<1e-12 and errors['force_samples/primary_force_N']<1e-18))
 write(BASE/'REGRESSION.json',reg);assert all(x['PASS'] for x in reg)
 for w in [1,2,4]:
  tasks=[task(f'workers{w}_{i}',duration=.00025) for i in range(4)];t=time.perf_counter();rs=batch(tasks,w);wall=time.perf_counter()-t
  parallel.append(dict(workers=w,wall_s=wall,trajectory_wall_sum_s=sum(r['wall_s'] for r in rs),PASS=all(r['status']=='PASS' for r in rs)));records+=rs
 # Four-amplitude G1 and phase usefulness G2 preregistered controls, including OFF.
 tasks=[task('g1_off',amps=[0]*4,duration=.0005),task('g1_x_half',amps=[0,0,500,500],duration=.0005),task('g1_y',amps=[700,700,0,0],phase=[0,np.pi/2,0,0],duration=.0005),task('g1_four',amps=[700]*4,phase=[0,np.pi/2,0,np.pi],duration=.0005),task('g2_phase',phase=[0,0,0,7*np.pi/6],duration=.0005)]
 records+=batch(tasks,2)
 c=context();x=c['x0'];s=make_schedule([300,500,700,800],[0,.8,1.7,2.8],.001);a=AuditEngine(x,c['physics'],c['response'],c['basis'],c['cavity'],policy=policy());a.control(s.segments[0].event());f=a.fields(x)
 t=make_schedule([300,500,700,800],[v+.731 for v in [0,.8,1.7,2.8]],.001);a.control(t.segments[0].event());g=a.fields(x)
 gauge={k:float(np.max(abs(f[k]-g[k]))) for k in ['primary_force_N','secondary_force_N','total_force_N']}
 result=dict(regression=reg,parallel=parallel,records=records,gauge_force_difference_N=gauge,workers_selected=min(parallel,key=lambda r:r['wall_s'])['workers'],wall_s=time.perf_counter()-start,versions={n:md.version(n) for n in ['numpy','scipy','h5py','warp-lang','ribs','pymoo','scikit-learn']},identity=identity())
 write(BASE/'PILOT.json',result)
if __name__=='__main__':main()

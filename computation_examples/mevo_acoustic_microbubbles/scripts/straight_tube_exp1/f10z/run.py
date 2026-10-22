"""Reuse F10S Campaign.run verbatim, injecting only paths, states and schedules."""
import os,time,traceback,threading,subprocess,resource,cProfile,pstats
import numpy as np,h5py
from scripts.straight_tube_exp1.f10s import run as inherited
from scripts.straight_tube_exp1.f10s.common import context
from scripts.straight_tube_exp1.f10s.admission import layered
from scripts.straight_tube_exp1.f10p.analyze import comparison
from mevo_mb.control.scheduled_pzt import Schedule
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .common import *
from .state import series

class Campaign(inherited.Campaign):
 def __init__(self,out):
  self.started=time.perf_counter();self.out=out;self.config=read(CONFIG);self.p=policy();self.ident=identity();self.records=[]
  verify_inputs();c=context();self.x0=c['x0'];self.physics=c['physics'];self.response=c['response'];self.cavity=c['cavity'];self.bases={'coarse':c['basis'],'fine':CachedWaterBasis(OLD/'cache/fine.h5')}
  with h5py.File(BASE/'INITIAL_STATES.h5') as h:self.states={k:h[k][:] for k in h}
  self.prior_s=sum(read(p)['wall_s'] for p in BASE.glob('*/CAMPAIGN.json'));self.prior_count=sum(len(read(p).get('runs',[])) for p in BASE.glob('*/CAMPAIGN.json'))
  self.route_budget=self.config['remote_budget_s']-self.prior_s;self.schedules={k:Schedule.from_record(v) for k,v in read(BASE/'RUN_SCHEDULES.json').items()}
  # The inherited orchestrator writes only through self.out; these two immutable
  # provenance files belong to F10Z. No historical source is edited or duplicated.
  inherited.BASE=BASE
 def trial(self,initial,action,duration=.002,variant='DT1'):
  name=f'{initial}__{action}__{round(duration*1e6)}us__{variant}'
  old=list(BASE.glob(f'*/runs/{name}/RESULT.json'))
  if old:return read(old[0])
  result=self.run(name,action,duration=duration,dt=.5e-6 if variant=='DT05' else 1e-6,cache='fine' if variant=='FINE' else 'coarse',initial=initial)
  return result
 def triple(self,initial,action,duration=.002):
  results=[self.trial(initial,action,duration,v) for v in ['DT1','DT05','FINE']]
  return assess(initial,action,duration)

def locate(initial,action,duration,variant):
 return next(BASE.glob(f'*/runs/{initial}__{action}__{round(duration*1e6)}us__{variant}/RESULT.json'),None)
def assess(initial,action,duration):
 files=[locate(initial,action,duration,v) for v in ['DT1','DT05','FINE']]
 if any(p is None for p in files):return dict(admitted=False,reason='MISSING_REFINEMENT',initial=initial,action=action,duration_s=duration)
 records=[read(p) for p in files]
 if any(r['status']!='PASS' for r in records):return dict(admitted=False,reason='INCOMPLETE_TRAJECTORY',initial=initial,action=action,duration_s=duration)
 comparisons=[]
 for p in files[1:]:
  s,_=comparison(files[0].parent,p.parent,'F10Z unchanged-physics refinement')
  with h5py.File(files[0].parent/'trajectory.h5') as h:base_n=len(h['time_s'])
  s['complete_window']=s['common_frames']==base_n and abs(s['common_end_s']-duration)<1e-13;comparisons.append(s)
 top=[]
 for p in files:
  with h5py.File(p.parent/'trajectory.h5') as h:t=h['time_s'][:];x=h['positions_m'][:]
  top.append(dict(initial_network_persists_entire_window=bool(np.all(series(t,x)['components']==1))))
 gates=layered(records[0],*comparisons,top,policy(),duration)
 numerical=all(gates[k]=='PASS' for k in ['BASIC_NUMERICS_PASS','GROUP_SHAPE_NUMERICALLY_ADMITTED','PARTICLE_ID_NUMERICALLY_ADMITTED']) and all(s['COM_m']<=policy()['group_COM_tolerance_m'] for s in comparisons)
 # Initial single-component survival is intentionally NOT an admission requirement
 # for a dispersed state. Its graph event agreement is assessed independently.
 out=dict(admitted=bool(numerical),initial=initial,action=action,duration_s=duration,comparisons=comparisons,legacy_layered=gates,topology_scope='Original single-network-survival gate reported; dispersed event agreement analyzed separately, not relabelled PASS')
 write(BASE/'checks'/f'{initial}__{action}__{round(duration*1e6)}us.json',out);return out

class Monitor:
 def __init__(self):self.done=threading.Event();self.rows=[];self.thread=threading.Thread(target=self.loop,daemon=True)
 def loop(self):
  start=time.perf_counter();previous=start;previous_cpu=time.process_time()
  while not self.done.wait(1):
   now=time.perf_counter();cpu=time.process_time()
   rss=int(Path('/proc/self/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')/2**20
   row=dict(elapsed_s=now-start,rss_MiB=rss,cpu_percent=100*(cpu-previous_cpu)/(now-previous))
   previous=now;previous_cpu=cpu
   try:
    s=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True,timeout=2);row['gpu_percent'],row['gpu_MiB']=[float(v) for v in s.strip().split(',')]
   except Exception as e:row['gpu_query_error']=repr(e)
   self.rows.append(row)

def main():
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['A','B','C'],required=True);ap.add_argument('--label',required=True);a=ap.parse_args()
 out=BASE/a.label;out.mkdir(parents=True,exist_ok=False);(BASE/'checks').mkdir(exist_ok=True);start=time.perf_counter();c=None;mon=Monitor();mon.thread.start();profile=cProfile.Profile();profile.enable();status='EXCEPTION'
 try:
  c=Campaign(out)
  if a.stage=='A':
   c.checks()
   from .diagnostics import run as diagnose
   diagnose(c)
   for initial in c.config['primary_states']:
    for action in ['OFF','PAIR']:c.trial(initial,action,.0005)
    for action in ['OFF','XM','XP','PAIR']:
     screen=locate(initial,'PAIR' if action=='PAIR' else 'OFF',.0005,'DT1')
     if read(screen)['status']=='PASS':c.trial(initial,action)
   # Original complete 3ms result is independently replayed, no concatenation.
   c.triple('original','PAIR',.003)
  elif a.stage=='B':
   for initial in c.config['primary_states']:c.triple(initial,'PAIR')
   for action in ['OFF','XM','XP']:c.triple('medium_train',action)
   for initial in c.config['heldout_states']:
    for action in ['OFF','PAIR']:c.triple(initial,action)
   for action in ['AMP-5','AMP+5','PHASE-3','PHASE+3','PHASE-6','PHASE+6']:c.triple('medium_train',action)
   for action in ['YM','YP']:c.triple('medium_train',action)
   c.triple('medium_train','PAIR_OFF')
  else:
   previous=.002;decisions=[]
   for duration in c.config['extension_s']:
    gate=assess('medium_train','PAIR',previous);p=locate('medium_train','PAIR',previous,'DT1')
    if p:
     with h5py.File(p.parent/'trajectory.h5') as h:m=series(h['time_s'][:],h['positions_m'][:])
     cost=read(p)['wall_s']*duration/previous*4.5;cost_ok=cost<180 and cost+c.prior_s+time.perf_counter()-c.started<c.config['remote_budget_s'];shape_cost=abs(m['rg_m'][-1]/m['rg_m'][0]-1)
    else:cost=float('inf');cost_ok=False;shape_cost=float('inf')
    eligible=gate['admitted'] and cost_ok and shape_cost<=.05
    decision=dict(previous_s=previous,next_s=duration,previous_gate=gate,projected_triple_wall_s=cost,relative_Rg_change=shape_cost,eligible=bool(eligible),reason='Numerical/refinement, cost and <=5% deformation scope; no validated bubble fusion physics')
    decisions.append(decision);write(out/'EXTENSION_DECISIONS.json',decisions)
    if not eligible:break
    g=c.triple('medium_train','PAIR',duration)
    if not g['admitted']:break
    previous=duration
  status='COMPLETED'
 except Exception:
  (out/'EXCEPTION.txt').write_text(traceback.format_exc());raise
 finally:
  profile.disable();profile.dump_stats(str(out/'profile.pstats'))
  with (out/'profile.txt').open('w') as f:pstats.Stats(profile,stream=f).sort_stats('cumtime').print_stats(65)
  mon.done.set();mon.thread.join(4);write(out/'resources.json',mon.rows)
  write(out/'CAMPAIGN.json',dict(stage=a.stage,status=status,wall_s=time.perf_counter()-start,runs=c.records if c else [],identity=identity(),CPU_user_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime,CPU_system_s=resource.getrusage(resource.RUSAGE_SELF).ru_stime))
if __name__=='__main__':main()

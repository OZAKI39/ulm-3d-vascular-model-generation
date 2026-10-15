from checks import *
from typed_model import TypedLPS
class Integrator:
 def __init__(self,ns,journal):
  self.ns=ns;self.j=journal;self.ref=ns.ref;self.model=TypedLPS(ns.ref);self.queue=[];self.trials=[];self.rollbacks=[];self.serial=0;self.accepts=0;self.rejects=0;self.base_intervals=0;self.native_wall=0.;self.restore_wall=0.
 def next(self,p):
  if self.queue:h,depth,start,end=self.queue.pop(0)
  else:h,depth,start,end=min(BASE,END-p['time']),0,p['time'],min(p['time']+BASE,END)
  while True:
   self.serial+=1;xp=predict_x(p,h,self.ref);_,dp=self.model.update(xp,p['kappa']);forecast=float(np.max(dp-p['d']));self.j.trial_start(self.serial,p,h,depth,None)
   assert np.isfinite(xp).all(),'PREDICTED_NONFINITE'
   checkpoint=forecast>=CHECKPOINT_TRIGGER or depth>0
   if checkpoint:
    cp=self.j.checkpoint(self.ns,p,'damage_trial');same=compare(p,self.ns.state());assert same['status']=='PASS',same;before_top=topology(self.ns)
   if h!=self.ns.dt:self.ns.l.command(f'timestep {h:.17g}');self.ns.dt=h
   ts=time.perf_counter();self.ns.l.command('run 1 pre no post no');s=self.ns.state();self.native_wall+=time.perf_counter()-ts
   delta=s['d']-p['d'];maxdd=float(delta.max());bad=np.flatnonzero(delta>.1);err=float(abs(s['x']-xp).max());record=dict(trial=self.serial,from_step=p['step'],candidate_step=s['step'],start_s=p['time'],end_s=s['time'],h_s=h,depth=depth,max_dd=maxdd,forecast_dd=forecast,predictor_x_error_m=err,decision='REJECT' if len(bad) else 'ACCEPT_DAMAGE_GATE',violating_ids=bad+1);self.trials.append(record);self.j.line(self.j.trials,dict(kind='DAMAGE_CHECK',**record))
   assert err<2e-18 and abs(s['time']-p['time']-h)<1e-15,(err,s['time'],h)
   assert s['step']==p['step']+1 and s['scalars'][4]==p['scalars'][4]+1 and s['scalars'][5]==s['step']
   if len(bad):
    self.rejects+=1;save_state(O/'full/rejected'/f'trial_{self.serial:06d}.npz',s);assert checkpoint,'Uncheckpointed damage violation'
    save_state(O/'full/rejected'/f'before_{self.serial:06d}.npz',p);ts=time.perf_counter();query=self.ns.ext.query;self.ns.close();self.ns=Native(f'full/restore_{self.rejects:03d}',restart=cp,restart_state=p,dt=h,query=query);restored=self.ns.dynamics();after_top=topology(self.ns);check=compare(p,restored);check.update(trial=self.serial,native_topology_exact=np.array_equal(before_top[0],after_top[0]),special_exclusions_exact=before_top[1]==after_top[1]);self.rollbacks.append(check);self.j.event('ROLLBACK',**check);self.restore_wall+=time.perf_counter()-ts
    assert check['status']=='PASS' and check['native_topology_exact'] and check['special_exclusions_exact'],check
    assert h>MIN_H*(1+1e-12),'MINIMUM_EVENT_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1;self.queue.insert(0,(h,depth,start,end));p=restored;continue
   self.accepts+=1
   if not self.queue:assert abs(s['time']-end)<1e-15;self.base_intervals+=1
   return s,h
 def save(self):dump(O/'full/TRIAL_ACCEPT_REJECT_AUDIT.json',dict(accepted_damage_gate=self.accepts,rejected=self.rejects,base_intervals=self.base_intervals,trials=self.trials,rollbacks=self.rollbacks,native_candidate_wall_s=self.native_wall,restore_wall_s=self.restore_wall))

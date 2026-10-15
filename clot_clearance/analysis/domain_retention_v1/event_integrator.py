"""Native VV candidates, native restart rollback, and dyadic physical substeps.

The predictor only selects checkpoints; actual native candidate damage decides.
Physical histories are never edited by this driver.
"""
from native import *
import time

def boundary(t):
 ql=min(t/.0003,1.);qu=0. if t<.00035 else min((t-.00035)/.0003,1.)
 return (30e-6*(np.cos(np.pi*qu)-np.cos(np.pi*ql)),
         60e-6*np.pi/(2*.0003)*(np.sin(np.pi*ql)-np.sin(np.pi*qu)))

def topology(ns):
 nb,raw=ns.l.gather_bonds();b=np.ctypeslib.as_array(raw,shape=(3*nb,)).reshape(-1,3)[:,1:].copy();b.sort(axis=1)
 codes=np.sort(b[:,0]*ns.N+b[:,1])
 fn=ns.lib.lps_event_special;fn.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_double)];fn.restype=ctypes.c_int
 count=fn(ns.l.lmp,None);out=np.empty(count);assert fn(ns.l.lmp,ns.ptr(out))==count
 exclusions=[];at=0
 while at<count:
  atom_id=int(out[at]);counts=out[at+1:at+4].astype(int);at+=4
  entries=out[at:at+counts[2]].astype(int);at+=counts[2]
  exclusions.append((atom_id,counts.tolist(),[sorted(entries[:counts[0]].tolist()),sorted(entries[counts[0]:counts[1]].tolist()),sorted(entries[counts[1]:counts[2]].tolist())]))
 exclusions.sort()
 return codes,exclusions

def restore_comparison(before,after):
 err={k:float(np.max(np.abs(np.asarray(before[k])-np.asarray(after[k])))) for k in ['x','velocity','force','material_force','kappa','d','time','step','native_bonds','Uc','K','native_U']}
 err['physical_scalars']=float(abs(before['scalars'][:6]-after['scalars'][:6]).max())
 passed=all(err[k]==0 for k in ['x','velocity','kappa','d','time','step','native_bonds'])
 passed &= max(err['force'],err['material_force'])<=2e-17 and max(err['Uc'],err['K'],err['native_U'],err['physical_scalars'])<=1e-20
 return dict(pass_=bool(passed),absolute_errors=err,diagnostic_evaluation_counters='Reinitialized by restart, not physical history; global attempts recorded separately')

class EventIntegrator:
 def __init__(self,ns,case,ids,journal=None):
  self.journal=journal
  self.ns=ns;self.case=case;self.folder=O/case;self.ids=ids;self.ref=ns.ref;self.model=SofteningLPS(self.ref)
  self.fi=ids['free_ids']-1;self.queue=[];self.trials=[];self.rollbacks=[];self.checkpoints=0;self.rejects=0;self.accepts=0;self.checkpoint_s=0.;self.restore_s=0.;self.predict_s=0.;self.native_s=0.;self.intervals=[];self.trial_serial=0
  (self.folder/'transactions').mkdir(exist_ok=True)
 def predict(self,p,h):
  start=time.perf_counter();x=p['x'].copy();fi=self.fi
  # Native nve applies the first half kick before drift.
  x[fi]+=h*(p['velocity'][fi]+h/(2*RHO*V)*p['force'][fi])
  u,_=boundary(p['time']+h)
  for key,sign in [('left_ids',-1),('right_ids',1)]:
   ii=self.ids[key]-1;x[ii]=self.ref['X'][ii];x[ii,0]+=sign*u
  _,d=self.model.update(x,p['kappa'],s_complete=SC)
  self.predict_s+=time.perf_counter()-start
  return x,float(np.max(d-p['d']))
 def next(self,p):
  if not self.queue:
   h=BASE;depth=0;interval_start=p['time'];interval_end=p['time']+BASE
  else:h,depth,interval_start,interval_end=self.queue.pop(0)
  while True:
   self.trial_serial+=1;serial=self.trial_serial
   xp,forecast=self.predict(p,h)
   clearance=float(np.min(np.minimum(xp+BOX,BOX-xp)))
   if self.journal:self.journal.trial_start(serial,p,h,depth,clearance)
   if clearance<=0:
    if self.journal:self.journal.checkpoint(self.ns,p,'predicted_domain_exit')
    raise RuntimeError('FROZEN_DOMAIN_WOULD_LOSE_ATOM: predicted candidate outside box')
   if clearance<DX and self.journal:self.journal.checkpoint(self.ns,p,'domain_proximity')
   checkpoint=forecast>=CHECKPOINT_TRIGGER or depth>0
   path=self.folder/'transactions'/f'accepted_before_trial_{serial:05d}.restart'
   if checkpoint:
    start=time.perf_counter();self.ns.l.command(f'write_restart {path}');self.checkpoint_s+=time.perf_counter()-start;self.checkpoints+=1
    check=restore_comparison(p,self.ns.state());assert check['pass_'],('checkpoint mutated state',check)
    tc,ec=topology(self.ns)
    if self.journal:self.journal.register_checkpoint(path,p,'damage_trial',self.ns)
   if self.ns.dt!=h:self.ns.l.command(f'timestep {h:.17g}');self.ns.dt=h
   start=time.perf_counter();self.ns.l.command('run 1 pre no post no');s=self.ns.state();self.native_s+=time.perf_counter()-start
   delta=s['d']-p['d'];maxdd=float(delta.max());bad=np.flatnonzero(delta>.10)
   pred_error=float(abs(s['x']-xp).max())
   record=dict(trial=serial,from_step=p['step'],candidate_step=s['step'],start_s=p['time'],end_s=s['time'],h_s=h,depth=depth,base_interval_start_s=interval_start,base_interval_end_s=interval_end,forecast_max_dd=forecast,max_dd=maxdd,max_bond_id=int(np.argmax(delta)+1),bond1813_kappa_before=p['kappa'][1812],bond1813_kappa_after=s['kappa'][1812],bond1813_d_before=p['d'][1812],bond1813_d_after=s['d'][1812],checkpoint=checkpoint,predictor_x_error_m=pred_error,decision='REJECT' if len(bad) else 'ACCEPT',violating_ids=bad+1)
   self.trials.append(record)
   if self.journal:self.journal.trial_result(record)
   assert pred_error<2e-18,('predictor coordinate mismatch',record)
   assert abs(s['time']-p['time']-h)<1e-15,('cumulative time mismatch',record)
   # Native virgin history has token=-1, commits=0. First physical token is
   # ntimestep=1, not token_old+1; thereafter the step identity remains exact.
   assert s['scalars'][4]==p['scalars'][4]+1 and s['scalars'][5]==s['step'] and s['step']==p['step']+1
   assert np.array_equal(self.ns.field(9),np.arange(1,self.ns.NB+1))
   if len(bad):
    self.rejects+=1
    np.savez_compressed(self.folder/'transactions'/f'rejected_trial_{serial:05d}.npz',**s)
    np.savez_compressed(self.folder/'transactions'/f'before_trial_{serial:05d}.npz',**p)
    assert checkpoint, 'Unsafe uncheckpointed candidate: stop, never accept'
    start=time.perf_counter();self.ns.close()
    self.ns=Native(f'{self.case}/restore_{self.rejects:03d}',restart=path,dt=h)
    self.ns.dynamics(self.ids['left_ids'],self.ids['right_ids'],self.ids['free_ids'])
    restored=self.ns.state();check=restore_comparison(p,restored);tc2,ec2=topology(self.ns)
    check.update(trial=serial,topology_exact=np.array_equal(tc,tc2),contact_exclusions_exact=ec==ec2,checkpoint=str(path.relative_to(O)),checkpoint_sha256=sha(path),driver_ledger_untouched=True)
    check['pass_'] &= check['topology_exact'] and check['contact_exclusions_exact']
    self.rollbacks.append(check);
    if self.journal:self.journal.event('rollback',check)
    self.restore_s+=time.perf_counter()-start
    np.savez_compressed(self.folder/'transactions'/f'restored_trial_{serial:05d}.npz',**restored)
    assert check['pass_'],('Incomplete native rollback',check)
    assert h>MIN_H*(1+1e-12),'MINIMUM_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1
    self.queue.insert(0,(h,depth,interval_start,interval_end))
    # Use the actual restored force, never a rejected candidate force.
    p=restored
    continue
   self.accepts+=1
   if depth:
    np.savez_compressed(self.folder/'transactions'/f'accepted_trial_{serial:05d}.npz',**s)
   if not self.queue:
    covered=s['time']-interval_start
    assert abs(covered-BASE)<1e-15,('Incomplete physical base interval',covered)
    self.intervals.append(dict(start_s=interval_start,end_s=s['time'],covered_s=covered))
   return s,h
 def save(self):
  dump(self.folder/'TRIAL_ACCEPT_REJECT_AUDIT.json',dict(trials=self.trials,accepted=self.accepts,rejected=self.rejects,checkpoints=self.checkpoints,base_intervals=self.intervals,all_intervals_exact_within_1e_15_s=all(abs(z['covered_s']-BASE)<1e-15 for z in self.intervals),checkpoint_wall_s=self.checkpoint_s,restore_wall_s=self.restore_s,predict_wall_s=self.predict_s,native_candidate_wall_s=self.native_s))
  dump(self.folder/'ROLLBACK_STATE_CONSISTENCY.json',dict(status='PASS' if all(z['pass_'] for z in self.rollbacks) else 'FAIL',records=self.rollbacks))

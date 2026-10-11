"""Frozen F7 replay M0 dt/2, with read-only accepted-step GPU diagnostics."""
import sys,json,hashlib,math
from pathlib import Path
import numpy as np
SNAP=Path('/workspace/mevo_acoustic_microbubbles/f7/.execution/9c198a0dd8589d6c6b5414ca529ab0a861cfe74e')
sys.path[:0]=[str(SNAP/'src'),str(SNAP)]
import warp as wp
from mevo_mb.fast_sim.f6_engine import F6Simulator
from mevo_mb.fast_sim import kernels as k
from scripts.fast_pzt_f7 import runner,common

@wp.kernel
def audit(x:wp.array(dtype=wp.vec3d),pred:wp.array(dtype=wp.vec3d),end:wp.array(dtype=wp.vec3d),counter:wp.array(dtype=int),stats:wp.array2d(dtype=wp.float64),n:int,R:wp.float64,dt:wp.float64):
 i=wp.tid();s=counter[0];v=(pred[i]-x[i])/dt
 wp.atomic_max(stats,s,0,wp.length(end[i]-x[i])/R)
 wp.atomic_max(stats,s,1,wp.length(pred[i]-x[i])/R)
 wp.atomic_max(stats,s,2,wp.length(end[i]-pred[i])/R)
 mingap=wp.float64(1.e9);closing=wp.float64(1.e9);approach=wp.float64(0.0);ratio=wp.float64(0.0)
 for j in range(n):
  if j!=i:
   d=wp.length(x[j]-x[i]);gap=d-wp.float64(2.0)*R
   endgap=wp.length(end[j]-end[i])-wp.float64(2.0)*R;mingap=wp.min(mingap,endgap)
   vj=(pred[j]-x[j])/dt
   closing_v=-wp.dot(vj-v,(x[j]-x[i])/d)
   approach=wp.max(approach,closing_v)
   if gap>wp.float64(0.0) and closing_v>wp.float64(0.0):
    tc=gap/closing_v;closing=wp.min(closing,tc);ratio=wp.max(ratio,dt/tc)
 wp.atomic_min(stats,s,3,mingap/R);wp.atomic_max(stats,s,4,approach)
 wp.atomic_min(stats,s,5,closing);wp.atomic_max(stats,s,6,ratio)
 if mingap<=wp.float64(0.02)*R:wp.atomic_add(stats,s,7,wp.float64(1.0)/wp.float64(n))
 if mingap<-wp.float64(0.01)*R:wp.atomic_add(stats,s,8,wp.float64(1.0))
@wp.kernel
def advance_counter(counter:wp.array(dtype=int)):counter[0]=counter[0]+1

class AuditSimulator(F6Simulator):
 instance=None
 def __init__(self,*a,**kw):
  super().__init__(*a,**kw);AuditSimulator.instance=self
  self.audit_pred=wp.zeros(self.total,dtype=wp.vec3d,device=self.device)
  self.audit_counter=wp.zeros(1,dtype=int,device=self.device)
  initial=np.zeros((200010,9));initial[:,3]=1e9;initial[:,5]=1e9
  self.audit_initial=initial;self.audit_stats=wp.array(initial,dtype=wp.float64,device=self.device)
  wp.load_module(module=sys.modules[__name__],device=self.device)
 def reset(self):
  super().reset()
  if hasattr(self,'audit_counter'):self.audit_counter.zero_();self.audit_stats.assign(self.audit_initial)
 def launch(self,kernel,dim,args):
  if kernel is k.accept and hasattr(self,'audit_counter'):
   wp.launch(audit,dim=self.total,inputs=[self.x,self.audit_pred,args[0],self.audit_counter,self.audit_stats,self.n,self.p.radius_m,self._current_dt],device=self.device)
  result=super().launch(kernel,dim,args)
  if kernel is k.predict and hasattr(self,'audit_counter'):wp.copy(self.audit_pred,args[2])
  return result
 def step(self,dt):
  super().step(dt)
  if hasattr(self,'audit_counter'):wp.launch(advance_counter,dim=1,inputs=[self.audit_counter],device=self.device)

if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
 baseline=Path('/workspace/mevo_acoustic_microbubbles/f7/results/fast_pzt/f7/replay_v1/M0')
 config=json.loads((baseline/'CONFIG.json').read_text());summary=json.loads((baseline/'SUMMARY.json').read_text())
 protocol=json.loads((SNAP/'configs/fast_pzt/f7/protocol.json').read_text())
 plan=dict(baseline=str(baseline),source_snapshot=str(SNAP),reason='Contact-active 32-particle mature cluster, baseline correction/free path=0.396 and accepted overlap/R=0.00759; monitor persistence/morphology across prescribed controls, not de-novo formation',window_s=[0,1],same_initial_ids=True,initial_seed='Explicit archived initial coordinates; no new randomness',change='Only nominal dt10->5us; no adaptive displacement cap in this actual F6 fixed_event_aligned backend',extra_kernel='Read-only accepted proposal before original commit; original kernel sequence/order unchanged; diagnostic arrays never feed forces',criteria=dict(cluster_fraction_absolute=.05,timeweighted_aggregation_absolute=.05,accepted_penetration_over_R=.01,contact_pass_correction_over_R=.1),hard_contact_label_tolerance_over_R=.001,hard_contact_label_is_not_solver_stop_tolerance=True,integrator='Overdamped forward Euler prediction plus12 Jacobi PBD passes; no RK intermediate stages',baseline_missing='Actual accepted displacement and positive-gap closing rates not logged at every baseline step; baseline quality contains prediction/pass-correction/accepted overlap extrema; saved states are0.5ms',new_computation_count=1,maximum_extra_displacement_cap_runs=0,physics='Uncalibrated exploratory regularized pair near-field; no experimental cluster labels')
 common.write(out/'test_plan.json',plan)
 common.write(out/'SOURCE_IDENTITY.json',dict(identity=common.identity(),actual_f6_module=sys.modules[F6Simulator.__module__].__file__,baseline_config_sha256=hashlib.sha256((baseline/'CONFIG.json').read_bytes()).hexdigest(),baseline_trajectory_sha256=hashlib.sha256((baseline/'trajectory.h5').read_bytes()).hexdigest()))
 runner.F6Simulator=AuditSimulator
 result=runner.run_case(out/'dt_half',config['environment'],config['events'],protocol,common.objects(),dt=config['dt_s']/2,mobility=config['mobility'],gradient=config['gradient'],duration=1.,adaptive=config['adaptive'],feedback=config['feedback'])
 sim=AuditSimulator.instance;count=int(sim.audit_counter.numpy()[0]);stats=sim.audit_stats.numpy()[:count]
 np.savez_compressed(out/'dt_half/accepted_step_audit.npz',values=stats,columns=['actual_displacement_over_R','start_prediction_over_R','endpoint_net_projection_over_R','minimum_gap_over_R','max_relative_approach_m_s','min_positive_gap_closing_s','max_dt_over_positive_gap_closing','nearcontact_particle_fraction','nodes_beyond_registered_overlap'])
 common.write(out/'INSTRUMENTATION.json',dict(accepted_kernel_calls=count,reported_accepted_steps=result['actual_substep_count'],maximums=stats[:,[0,1,2,4,6,7,8]].max(0).tolist(),minimum_gap_over_R=float(stats[:,3].min()),minimum_positive_closing_s=float(stats[:,5].min()),extra_physics_changes=False))

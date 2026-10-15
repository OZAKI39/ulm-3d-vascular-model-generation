"""Minimum deterministic replay to obtain a genuine pre-event native restart."""
from event_integrator import *
import resource
def main():
 freeze=json.loads((O/'EVENT_INTEGRATOR_FREEZE.json').read_text())
 for n,h in freeze['driver_sha256'].items():assert sha(A/n)==h
 with (O/'replay/STARTED.json').open('x') as f:json.dump(dict(purpose='native accepted 640.25 us checkpoint',pid=os.getpid()),f)
 ids=np.load(O/'provenance/341_boundary_ids.npz');ns=Native('replay/native',dt=BASE);ns.dynamics(ids['left_ids'],ids['right_ids'],ids['free_ids'])
 start=time.perf_counter();maxdd=0.;prev=ns.state()
 for step in range(1,2562):
  ns.l.command('run 1 pre no post no');s=ns.state();dd=float(np.max(s['d']-prev['d']));maxdd=max(maxdd,dd)
  assert dd<=.10 and s['scalars'][4]==step and s['scalars'][5]==step
  prev=s
 old=np.load(PREV/'keyframes/step_002561.npz');check=restore_comparison(old,s)
 assert check['pass_'],check
 ns.l.command(f'write_restart {O/"replay/accepted_640p25.restart"}')
 np.savez_compressed(O/'replay/accepted_640p25.npz',**s)
 dump(O/'replay/DETERMINISTIC_REPLAY_AUDIT.json',dict(status='PASS',comparison_to_previous_accepted_FP64=check,steps=2561,time_s=s['time'],max_dd=maxdd,wall_s=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,not_a_second_full_protocol=True,checkpoint_is_native=True))
 ns.close();print('REPLAY_PASS',flush=True)
if __name__=='__main__':main()

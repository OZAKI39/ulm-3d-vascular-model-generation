"""Auxiliary observer for the SAME running first CFD, no solver restart.
Corrects a read-only precision identity check and records native process exit.
"""
from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.measure_b import Measurements
from flow_solver_support.steady_monitor import SteadyStopMonitor
from flow_solver_support.transient_checks import checkpoint_audit
from flow_solver_support.flow_parser import parse_solver_log
from flow_solver_support.solver_checks import linear_gate,nonlinear_gate
import os,time,traceback
out=B/'baseline';run=out/'run';results=run/'8-procs';design=json.loads((B/'DESIGN_FREEZE.json').read_text());policy={k:design['criteria'][k] for k in ['save_interval_steps','steady_last_intervals','velocity_change_limit','flow_change_limit','mass_limit']};measure=Measurements(out/'SV_MESH/mesh_arrays.npz',design['flow']['Q_target_m3_s'],.005);monitor=SteadyStopMonitor(measure,policy);seen=set();stop=None;retries=[];started=json.loads((out/'FORMAL_ATTEMPT.json').read_text())['started_unix']
pids=[]
for d in Path('/proc').iterdir():
 if not d.name.isdigit():continue
 try:
  if (d/'comm').read_text().strip()=='mpirun' and os.readlink(d/'cwd')==str(run):pids.append(int(d.name))
 except OSError:pass
assert len(pids)==1;pidx=pids[0];exit_status=None
save(out/'MONITOR_ADAPTER_IDENTITY.json',dict(mpirun_pid=pidx,formal_CFD_attempt=1,solver_restart=False,reason='Native vtkPoints float32 but fields float64; exact bijective float32 roundtrip replaces overly strict geometric equality. Convergence thresholds unchanged.',reader_sha=sha(ROOT/'analysis/dual_track_v1/measure_b.py'),started_unix=time.time()))
def observe(done=False):
 global stop
 for p in sorted(results.glob('result_*.vtu')):
  step=int(p.stem.rsplit('_',1)[1]);cp=p.with_name('stFile_%03d.bin'%step)
  if step==0 or step%10 or step in seen or not cp.exists() or (not done and time.time()-max(cp.stat().st_mtime,p.stat().st_mtime)<2):continue
  try:checkpoint_audit(cp,step,.01,ranks=8);u,pressure=measure.read(p)
  except Exception as e:retries.append(dict(step=step,error=repr(e)));break
  seen.add(step);state=dict(measure.measure(u,pressure),step=step,time_s=step*.01)
  if monitor.observe(state,u) and stop is None:
   stop=dict(step=step,reason='Five consecutive steady intervals and finite/mass/wall gates, identical frozen criteria');(run/'STOP_SIM').write_text('0\n')
  save(out/'ADAPTER_PROGRESS.json',dict(states=monitor.states,intervals=monitor.intervals,stop=stop,elapsed_s=time.time()-started,read_retries=retries));print('observed',step,'mass',state['epsilon_mass'],'stop',stop,flush=True)
last=0
while True:
 try:
  raw=Path('/proc')/str(pidx)/'stat';text=raw.read_text();v=text[text.rfind(')')+2:].split()
  if v[0]=='Z':exit_status=os.waitstatus_to_exitcode(int(v[49]));save(out/'NATIVE_EXIT_STATUS.json',dict(mpirun_pid=pidx,raw_wait_status=int(v[49]),exit_code=exit_status,source='Own MPI launcher /proc/PID/stat zombie exit_code observed before parent reap',time_unix=time.time()));break
 except FileNotFoundError:break
 if time.monotonic()-last>3:observe();last=time.monotonic()
 time.sleep(.05)
observe(True)
# Native process is finished. Postprocess all original logs, no new flow solve.
log=(run/'solver.log').read_text(errors='replace');history=parse_solver_log(log,.01);save(out/'SOLVER_HISTORY.json',history);fail=[]
try:linear_gate(dict(exit_code=exit_status,history=history))
except Exception as e:fail.append(dict(gate='linear_or_exit',error=repr(e)))
try:nonlinear_gate(history,max_iterations=24)
except Exception as e:fail.append(dict(gate='nonlinear',error=repr(e)))
final=max(results.glob('result_*.vtu'),key=lambda p:int(p.stem.rsplit('_',1)[1]));u,p=measure.read(final);metrics=measure.measure(u,p);checks=dict(solver=exit_status==0 and not fail,steady=stop is not None,mass=max(metrics['epsilon_Q'],metrics['epsilon_mass'])<=1e-5,finite=metrics['velocity_finite'] and metrics['pressure_finite'],wall=metrics['wall_noslip_pass'],forward_outlets=all(v>0 for v in metrics['outlet_flows_m3_s'].values()))
save(out/'COORDINATE_IDENTITY_AUDIT.json',measure.identity);save(out/'CFD_MASS_BALANCE.json',dict(source='actual_vtu_surface_integration_on_frozen_double_input_geometry',**metrics));save(out/'CFD_RESULT.json',dict(status='PASS_WITH_SCOPE' if all(checks.values()) else 'FAIL',checks=checks,metrics=metrics,exit_code=exit_status,health_failures=fail,stop=stop,runtime_s=time.time()-started,snapshot=str(final.relative_to(ROOT)),snapshot_sha256=sha(final),mesh_sha256=sha(out/'SV_MESH/mesh-complete.mesh.vtu'),geometry_hash=design['geometry_sha256'],monitor_adapter=True,formal_CFD_attempts=1,solver_restart=False,scope='Steady SI water rigid V3 geometry, no streaming or FSI; mesh sensitivity separate'))

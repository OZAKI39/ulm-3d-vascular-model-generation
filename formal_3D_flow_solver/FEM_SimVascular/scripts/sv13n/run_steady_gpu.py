"""One GPU development run with the unchanged production steady monitor.

Native output every step permits STOP_SIM=0 without advancing to another output
boundary. Scientific monitor cadence remains exactly the production cadence.
No CPU comparison, profiling campaign, repeat timing, or signals to the solver.
"""
import json,math,re,subprocess,sys,time,hashlib
from pathlib import Path
from remote import ROOT,REMOTE,python,fetch
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import SteadyStopMonitor
from sv_validation.sv13n import production_monitor_policy_gate,checkpoint_one_rank
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv11 import linear_gate,nonlinear_gate
from sv_validation.provenance import sha256
from flow_parser import parse_solver_log
R=ROOT/'reports/sv1_3n';O=ROOT/'outputs/sv1_3n';S=ROOT/'scripts/sv13n'
name='REAL_VASCULAR_GPU';case=O/name;case.mkdir(exist_ok=True)
strategy=json.loads((ROOT/'configs/sv1_3n/convergence_strategy.json').read_text())
assert strategy['strategy']=='GPU_DEVELOPMENT_ONLY'
policy=json.loads((ROOT/strategy['production_policy_path']).read_text())
assert sha256(ROOT/strategy['production_policy_path'])==strategy['production_policy_sha256']
production_monitor_policy_gate(strategy['steady_policy'],policy)
assert json.loads((R/'svmp_gpu_smoke.json').read_text())['status']=='PASS'
ref=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],policy['Umean_m_s'])
monitor=SteadyStopMonitor(measure,policy);seen=set();stop=None;failure=None;start=time.monotonic()
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def remote(code):
 r=python(code,timeout=60);r.check_returncode();return json.loads(r.stdout)
def request_stop(kind,details):
 global stop
 if stop is not None:return
 result=remote(f'''from pathlib import Path
import json,time,os
c=Path({(REMOTE+'/outputs/'+name)!r});p=c/'STOP_SIM.pending';p.write_text('0\\n');os.replace(p,c/'STOP_SIM')
print(json.dumps(dict(request_unix_s=time.time(),mechanism='native STOP_SIM=0; finish only current timestep, write native restart and per-step VTU; no signals')))
''')
 stop=dict(kind=kind,detail=details,elapsed_s=time.monotonic()-start,**result)
 write('gpu_steady_stop_request',stop);print('SAFE_STOP '+json.dumps(stop),flush=True)
controller_log=ROOT/'logs/sv1_3n/steady_native_controller.log';assert not controller_log.exists()
with controller_log.open('x') as out:
 proc=subprocess.Popen([sys.executable,'-B',S/'invoke.py','flow_remote.py',name,'1','GPU','gpu13'],stdout=out,stderr=subprocess.STDOUT)
 while True:
  try:
   snapshot=remote(f'''from pathlib import Path
import json,time,re
b=Path({REMOTE!r});c=b/'outputs/{name}';now=time.time()
if c.exists():(c/'controller_heartbeat').write_text(str(now))
files=[]
for p in sorted((c/'1-procs').glob('result_*.vtu')):
 step=int(p.stem.rsplit('_',1)[1])
 if step%{policy['save_interval_steps']} or step in {sorted(seen)!r}:continue
 if now-p.stat().st_mtime<1:continue
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-200));closed=b'</VTKFile>' in f.read()
 cp=p.with_name('stFile_%03d.bin'%step)
 if closed and cp.exists():files.append(dict(step=step,path=str(p.relative_to(b)),checkpoint=str(cp.relative_to(b))))
log=b/'logs/{name}.log';last=None
if log.exists():
 with log.open('rb') as f:f.seek(max(0,log.stat().st_size-262144));tail=f.read().decode(errors='replace')
 rows=re.findall(r'^\\s*NS\\s+\\d+-.*$',tail,re.M);last=rows[-1] if rows else None
print(json.dumps(dict(files=files,last=last,execution_exists=(b/'reports/{name}_execution.json').exists())))
''')
   for f in snapshot['files']:
    path=O/Path(f['path']).relative_to('outputs');cp=O/Path(f['checkpoint']).relative_to('outputs')
    fetch(f['path'],path);fetch(f['checkpoint'],cp)
    fetch('logs/'+name+'.log',ROOT/'logs/sv1_3n'/('monitor_'+name+'.log'))
    h=parse_solver_log((ROOT/'logs/sv1_3n'/('monitor_'+name+'.log')).read_text(errors='replace'),policy['dt_s'])
    # A saved field is accepted only after the corresponding complete solver row.
    rows=[r for r in h['linear_solves'] if r['step']<=f['step']]
    assert rows and max(r['step'] for r in rows)==f['step']
    relevant=dict(h,linear_solves=rows,petsc_reasons=[r['petsc_reason'] for r in rows if r.get('petsc_reason')],unassigned_petsc_monitor=[],unassigned_petsc_reason=None)
    linear_gate(dict(history=relevant,exit_code=0,monitor_stop=None))
    nonlinear_gate(relevant,max_iterations=12,tolerance=1e-10,min_iterations=2)
    checkpoint=checkpoint_one_rank(cp,f['step'],policy['dt_s'])
    u,p=measure.read(path);state=measure.measure(u,p)
    state.update(step=f['step'],time_s=f['step']*policy['dt_s'],path=str(path.relative_to(ROOT)),sha256=sha256(path),checkpoint=checkpoint)
    assert state['velocity_finite'] and state['pressure_finite'] and state['wall_noslip_pass']
    assert all(math.isfinite(v) for v in [state['Q_in_m3_s'],*state['outlet_flows_m3_s'].values(),state['epsilon_mass']])
    eligible=monitor.observe(state,u,0,0);seen.add(f['step'])
    write('gpu_steady_history',dict(status='RUNNING',states=monitor.states,intervals=monitor.intervals,stop_eligible=eligible,policy_sha256=strategy['production_policy_sha256'],scientific_equivalence='DEFERRED',formal_check_cadence_steps=policy['save_interval_steps'],VTU_output_cadence_steps=1))
    print(f"GPU step {f['step']}: mass={state['epsilon_mass']:.6g}; steady={eligible}",flush=True)
    if eligible and stop is None:request_stop('STEADY',dict(first_full_steady_step=f['step'],checkpoint_sha256=checkpoint['sha256'],required_consecutive_intervals=policy['steady_last_intervals']))
   write('gpu_running',dict(status='RUNNING' if proc.poll() is None else 'EXITED',elapsed_s=time.monotonic()-start,last=snapshot['last'],saved_monitor_steps=sorted(seen),stop=stop,failure=failure))
  except Exception as exc:
   failure=failure or dict(reason='MONITOR_OR_FIELD_ERROR',detail=repr(exc))
   write('gpu_monitor_failure',failure)
   request_stop('FAILURE',failure)
  if proc.poll() is not None:break
  time.sleep(2)
 code=proc.wait()
subprocess.run([sys.executable,'-B',S/'fetch_case.py',name],check=True)
subprocess.run([sys.executable,'-B',S/'accept_case.py',name],check=True)
accept=json.loads((R/(name+'_acceptance.json')).read_text())
assert code==0 and failure is None and stop and stop['kind']=='STEADY','GPU_STEADY_NOT_REACHED'
final=accept['steps_completed'];vtu=case/f'1-procs/result_{final:03d}.vtu';cp=case/f'1-procs/stFile_{final:03d}.bin'
checkpoint=checkpoint_one_rank(cp,final,policy['dt_s'])
reload=json.loads(subprocess.check_output([sys.executable,'-B',S/'reload_field.py',vtu],text=True));assert reload['velocity_finite'] and reload['pressure_finite']
u,p=measure.read(vtu);state=measure.measure(u,p)
assert state['wall_noslip_pass'] and max(state['epsilon_Q'],state['epsilon_mass'])<=policy['mass_limit']
summary=dict(status='PASS',classification='GPU_STEADY_CANDIDATE',scientific_equivalence='DEFERRED',first_full_steady_step=stop['detail']['first_full_steady_step'],stop_step=final,physical_time_s=final*policy['dt_s'],required_consecutive_intervals=policy['steady_last_intervals'],last_interval=monitor.intervals[-1],stop_request=stop,safe_stop=True,reload=reload,checkpoint=checkpoint,measurement=state,VTU=str(vtu.relative_to(ROOT)),VTU_sha256=sha256(vtu),config_sha256=sha256(case/'solver.xml'),mesh_sha256=sha256(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz'),solver_sha256=accept['solver_sha256'],PETSc_library_sha256=accept['PETSc_library_sha256'],wall_time_s=accept['wall_time_s'],steps=final,linear_solves=accept['linear_solves'],total_KSP_iterations=sum(r['linear_iterations'] for r in accept['history']['linear_solves']),production='CPU_EARLY_STOP_PRODUCTION',production_changed=False,normal_exit=accept['exit_code']==0)
write('gpu_steady_candidate',summary)
history=json.loads((R/'gpu_steady_history.json').read_text());history.update(status='PASS',first_full_steady_step=summary['first_full_steady_step'],stop_step=final);write('gpu_steady_history',history)
print('GPU_STEADY_CANDIDATE complete; all further CFD/science/performance work deferred.',flush=True)

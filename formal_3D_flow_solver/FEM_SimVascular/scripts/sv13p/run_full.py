"""The only Stage P full steady run; original production monitor, low-frequency I/O."""
import json,math,re,subprocess,sys,time
from pathlib import Path
from remote import ROOT,REMOTE,python,fetch,upload
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import SteadyStopMonitor
from sv_validation.sv13n import production_monitor_policy_gate,checkpoint_one_rank
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv11 import linear_gate,nonlinear_gate
from sv_validation.provenance import sha256
from flow_parser import parse_solver_log
R=ROOT/'reports/sv1_3p';C=ROOT/'configs/sv1_3p';S=ROOT/'scripts/sv13p';O=ROOT/'outputs/sv1_3p'
name='REAL_VASCULAR_GPU_PC_WINNER';case=O/name
winner=json.loads((R/'winner.json').read_text());assert winner['full_required'];policyO=json.loads((C/'policy.json').read_text())
policy=json.loads((ROOT/policyO['production_policy_path']).read_text())
assert sha256(ROOT/policyO['production_policy_path'])==policyO['production_policy_sha256']
production_monitor_policy_gate(policyO['production_policy'],policy)
plan=dict(name=name,mode='full',candidate=winner['candidate'],MPI_ranks=1,GPUs=1,OMP_NUM_THREADS=1,start_step=0,end_step=policy['maximum_total_steps'],PETSC_OPTIONS=winner['PETSC_OPTIONS'],build_report=winner['build_report'])
p=C/'runplans'/(name+'.json');assert not p.exists();p.write_text(json.dumps(plan,indent=2)+'\n');upload(p,'configs/runplans/'+p.name)
ref=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],policy['Umean_m_s'])
monitor=SteadyStopMonitor(measure,policy);seen=set();stop=None;failure=None;start=time.monotonic()
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def remote(code):
 for attempt in range(3):
  try:
   r=python(code,timeout=25);r.check_returncode();return json.loads(r.stdout)
  except (subprocess.SubprocessError,json.JSONDecodeError):
   if attempt==2:raise
   time.sleep(attempt+1)
def stop_native(kind,detail):
 global stop
 if stop is not None:return
 observed=remote(f"from pathlib import Path\nimport json,time,os\nc=Path({(REMOTE+'/outputs/'+name)!r});p=c/'STOP_SIM.pending';p.write_text('0\\n');os.replace(p,c/'STOP_SIM')\nprint(json.dumps(dict(request_unix_s=time.time())))")
 stop=dict(kind=kind,detail=detail,elapsed_s=time.monotonic()-start,mechanism='Atomic native STOP_SIM=0; finish current timestep and force final VTU + full native checkpoint',**observed)
 write('steady_stop_request',stop);print('SAFE_STOP '+json.dumps(stop),flush=True)
controller=ROOT/'logs/sv1_3p/full_controller.log';assert not controller.exists()
with controller.open('x') as out:
 proc=subprocess.Popen([sys.executable,'-B',S/'invoke.py','flow_remote.py',name],stdout=out,stderr=subprocess.STDOUT)
 while True:
  try:
   snap=remote(f'''from pathlib import Path
import json,re,time
b=Path({REMOTE!r});c=b/'outputs/{name}';now=time.time()
if c.exists():(c/'controller_heartbeat').write_text(str(now))
files=[]
for p in sorted((c/'1-procs').glob('result_*.vtu')):
 step=int(p.stem.rsplit('_',1)[1])
 if step%10 or step in {sorted(seen)!r} or now-p.stat().st_mtime<1:continue
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-200));closed=b'</VTKFile>' in f.read()
 cp=p.with_name('stFile_%03d.bin'%step)
 if closed and cp.exists():files.append(dict(step=step,path=str(p.relative_to(b)),checkpoint=str(cp.relative_to(b))))
p=b/'logs/{name}.log';last=None
if p.exists():
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-262144));tail=f.read().decode(errors='replace')
 rows=re.findall(r'^\\s*NS\\s+\\d+-.*$',tail,re.M);last=rows[-1] if rows else None
print(json.dumps(dict(files=files,last=last)))
''')
   for f in snap['files']:
    vtu=O/Path(f['path']).relative_to('outputs');cp=O/Path(f['checkpoint']).relative_to('outputs')
    fetch(f['path'],vtu);fetch(f['checkpoint'],cp)
    logfile=ROOT/'logs/sv1_3p/monitor_full.log';fetch('logs/'+name+'.log',logfile)
    h=parse_solver_log(logfile.read_text(errors='replace'),policy['dt_s']);rows=[r for r in h['linear_solves'] if r['step']<=f['step']]
    assert rows and max(r['step'] for r in rows)==f['step']
    relevant=dict(h,linear_solves=rows,petsc_reasons=[r['petsc_reason'] for r in rows if r.get('petsc_reason')],unassigned_petsc_monitor=[],unassigned_petsc_reason=None)
    linear_gate(dict(history=relevant,exit_code=0,monitor_stop=None));nonlinear_gate(relevant,12,1e-10,2)
    checkpoint=checkpoint_one_rank(cp,f['step'],policy['dt_s']);u,p=measure.read(vtu);state=measure.measure(u,p)
    assert state['velocity_finite'] and state['pressure_finite'] and state['wall_noslip_pass']
    assert all(math.isfinite(v) for v in [state['Q_in_m3_s'],*state['outlet_flows_m3_s'].values(),state['epsilon_mass']])
    state.update(step=f['step'],time_s=f['step']*policy['dt_s'],path=str(vtu.relative_to(ROOT)),sha256=sha256(vtu),checkpoint=checkpoint)
    eligible=monitor.observe(state,u,0,0);seen.add(f['step'])
    write('steady_history',dict(status='RUNNING',states=monitor.states,intervals=monitor.intervals,stop_eligible=eligible,policy_sha256=policyO['production_policy_sha256'],scientific_equivalence='DEFERRED',VTU_cadence=10,formal_check_cadence=10))
    print(f"GPU PC WINNER step {f['step']}: mass={state['epsilon_mass']:.6g}, steady={eligible}",flush=True)
    if eligible and stop is None:stop_native('STEADY',dict(first_full_steady_step=f['step'],required_consecutive_intervals=5,checkpoint_sha256=checkpoint['sha256']))
   write('full_running',dict(status='RUNNING' if proc.poll() is None else 'EXITED',elapsed_s=time.monotonic()-start,last=snap['last'],saved_monitor_steps=sorted(seen),stop=stop,failure=failure))
  except Exception as exc:
   failure=failure or dict(reason='MONITOR_OR_FIELD_ERROR',detail=repr(exc));write('full_monitor_failure',failure);stop_native('FAILURE',failure)
  if proc.poll() is not None:break
  time.sleep(2)
 code=proc.wait()
subprocess.run([sys.executable,'-B',S/'fetch_case.py',name],check=True)
subprocess.run([sys.executable,'-B',S/'accept_case.py',name],check=True)
accepted=json.loads((R/(name+'_acceptance.json')).read_text())
assert code==0 and failure is None and stop and stop['kind']=='STEADY'
assert accepted['status']=='PASS' and max(accepted['measurement']['epsilon_Q'],accepted['measurement']['epsilon_mass'])<=policy['mass_limit']
baseline=json.loads((R/'reference_freeze.json').read_text())['GPU_baseline'];final=accepted['stop_step']
candidate=dict(status='PASS',classification='GPU_PC_WINNER_STEADY_CANDIDATE',scientific_equivalence='DEFERRED',first_full_steady_step=stop['detail']['first_full_steady_step'],stop_step=final,physical_time_s=final*policy['dt_s'],last_interval=monitor.intervals[-1],required_consecutive_intervals=5,safe_stop=True,normal_exit=True,reload=accepted['reload'],checkpoint=accepted['checkpoint'],measurement=accepted['measurement'],VTU=accepted['reload']['path'],VTU_sha256=accepted['reload']['sha256'],config_sha256=accepted['config_sha256'],solver_sha256=accepted['solver_sha256'],PETSc_library_sha256=accepted['PETSc_library_sha256'],mesh_sha256=sha256(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz'),winner=winner['candidate'],PETSC_OPTIONS=winner['PETSC_OPTIONS'],wall_time_s=accepted['wall_time_s'],steps=accepted['steps'],statistics=accepted['statistics'],VTU_count=accepted['VTU_count'],VTU_steps=accepted['VTU_steps'],sampled_device_memory_max_MiB=accepted['sampled_device_memory_max_MiB'],development_speedup=baseline['wall_time_s']/accepted['wall_time_s'],speedup_label='OBSERVATIONAL DEVELOPMENT SPEEDUP — SINGLE-RUN DEVELOPMENT COMPARISON',CPU_production='CPU_EARLY_STOP_PRODUCTION',production_changed=False)
write('winner_steady_candidate',candidate)
history=json.loads((R/'steady_history.json').read_text());history.update(status='PASS',first_full_steady_step=candidate['first_full_steady_step'],stop_step=final);write('steady_history',history)
print('Unique optimized full steady run complete; no further CFD.',flush=True)

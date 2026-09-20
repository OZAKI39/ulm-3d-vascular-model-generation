#!/usr/bin/env python3
"""Native CPU early stop. Never kills, suspends or edits the solver process."""
import json,os,re,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import *
from sv_validation.sv11 import ROW,parse_solver_log,linear_gate,nonlinear_gate,parse_boundary,write_csv
from sv_validation.sv12 import checkpoint_audit
from sv_validation.provenance import write_json,now
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.validation import integral_agreement

def stats(values):
    return {'count':len(values),'total':int(sum(values)),'mean':float(np.mean(values)),
            'median':float(np.median(values)),'P95':float(np.quantile(values,.95)),'max':int(max(values))}

check_reference();assert load('steady_stop_replay')['status']=='PASS'
policy=frozen_policy();ref=load('reference_freeze');opts=json.loads((CONFIG/'petsc_options.json').read_text())
assert opts==json.loads((ROOT/'configs/sv1_2/petsc_options.json').read_text())
build=load('petsc_build_manifest','sv1_1');case=OUTPUT/'cpu_early_stop'
assert sha256(case/'solver.xml')==sha256(ROOT/'configs/sv1_2/sv_flow.xml')
assert sha256(case/'4-procs/stFile_last.bin')==load('cpu_initial_state')['sha256']
log=LOG/'cpu_early_stop.log';resource=LOG/'cpu_early_stop_resources.txt'
assert not log.exists() and not (case/'STOP_SIM').exists(),'Never overwrite an experiment'
for path in (Path.home()/'.petscrc',case/'.petscrc',case/'petscrc'):assert not path.exists()
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],policy['Umean_m_s'])
monitor=SteadyStopMonitor(measure,policy);rows=[];offset=0;carry='';pending=[];failure=None;stop_request=None;stable={}
initial_step=load('cpu_initial_state')['records'][0]['step'];save=policy['save_interval_steps']

def completed_nonlinear_failures(step):
    final={r['step']:r for r in rows if r['step']<=step}
    return sum(not (r['nonlinear_iteration']>=2 and r['residual_finite'] and
        min(r['nonlinear_Ri_over_R0'],r['nonlinear_Ri_over_R1'])<=1e-10) for r in final.values())

def observe(path):
    u,p=measure.read(path);state=measure.measure(u,p);step=int(path.stem.rsplit('_',1)[1])
    state.update(step=step,time_s=step*policy['dt_s'],path=str(path.relative_to(ROOT)),sha256=sha256(path),
                 origin='inherited native seed' if step==initial_step else 'new native CPU run')
    if step>initial_step:require(any(r['step']==step for r in rows),'Saved field missing complete solver history')
    good=monitor.observe(state,u,sum(not r['linear_converged'] for r in rows),completed_nonlinear_failures(step))
    write_json(REPORT/'cpu_saved_states.json',{'states':monitor.states,'intervals':monitor.intervals,'stop_eligible':good})
    print(f"CPU saved {step}: mass={state['epsilon_mass']:.4g}; stop_eligible={good}",flush=True)
    return good

observe(case/f'4-procs/result_{initial_step:03d}.vtu')
seen={initial_step};last_status=0
env=dict(os.environ,PATH='/usr/bin:/bin',LD_LIBRARY_PATH=build['runtime_library_path'],
         PETSC_OPTIONS=opts['PETSC_OPTIONS'],OMP_NUM_THREADS=str(opts['OMP_NUM_THREADS']))
for key in ('PETSC_OPTIONS_YAML','PETSC_OPTIONS_TABLE','PETSC_OPTIONS_FILE'):env.pop(key,None)
cmd=['/usr/bin/mpiexec','--bind-to','none','-n',str(opts['mpi_ranks']),build['executable'],'solver.xml']
start=time.monotonic()
with log.open('x') as out, (OUTPUT/'qc/cpu_linear_history.jsonl').open('x') as stream:
    proc=subprocess.Popen(['/usr/bin/time','-v','-o',str(resource),*cmd],cwd=case,env=env,stdout=out,stderr=subprocess.STDOUT)
    def request_stop(kind,detail):
        global stop_request
        if stop_request:return
        latest=max([initial_step]+[r['step'] for r in rows]+list(seen))
        target=(latest//save+1)*save
        temp=case/'STOP_SIM.pending';temp.write_text(str(target)+'\n');temp.replace(case/'STOP_SIM')
        stop_request={'kind':kind,'detail':detail,'observed_step':latest,'requested_final_step':target,
                      'utc':now(),'elapsed_s':time.monotonic()-start,'mechanism':'native STOP_SIM on next save boundary; no signals'}
        write_json(REPORT/'cpu_stop_request.json',stop_request)
        print('Native safe stop requested: '+json.dumps(stop_request),flush=True)
    while True:
        alive=proc.poll() is None
        with log.open() as inp:inp.seek(offset);chunk=inp.read();offset=inp.tell()
        joined=carry+chunk
        if '\n' in joined:complete,carry=joined.rsplit('\n',1)
        else:complete='';carry=joined
        for line in complete.splitlines():
            pending.append(line)
            if re.search(r'DIVERGED_\w+|PETSC ERROR|ill-conditioned LHS|Resetting restart flag|(?<![a-z])(?:nan|inf)(?![a-z])',line,re.I):
                failure=failure or {'reason':'SOLVER_FAILURE','detail':line};request_stop('FAILURE',failure)
            if ROW.match(line):
                row=parse_solver_log('\n'.join(pending),policy['dt_s'])['linear_solves'][-1]
                row['linear_solve_index']=len(rows)+1;rows.append(row);pending=[]
                stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                if (len(rows)==1 and row['step']!=initial_step+1) or not row['linear_converged'] or not row.get('petsc_reason'):
                    failure=failure or {'reason':'INVALID_LINEAR_HISTORY','detail':row};request_stop('FAILURE',failure)
                if row['nonlinear_iteration']>=12 and min(row['nonlinear_Ri_over_R0'],row['nonlinear_Ri_over_R1'])>1e-10:
                    failure=failure or {'reason':'NONLINEAR_CONVERGENCE_FAIL','step':row['step']};request_stop('FAILURE',failure)
        for path in sorted((case/'4-procs').glob('result_*.vtu')):
            step=int(path.stem.rsplit('_',1)[1])
            if step in seen:continue
            s=path.stat();key=(s.st_size,s.st_mtime_ns);old=stable.get(step)
            stable[step]=(key,time.monotonic() if old is None or old[0]!=key else old[1])
            if not alive or (old and old[0]==key and time.monotonic()-old[1]>=2):
                try:
                    good=observe(path);seen.add(step)
                    cp=checkpoint_audit(case/f'4-procs/stFile_{step:03d}.bin',step,policy['dt_s'])
                    if good and failure is None:request_stop('STEADY',{'qualifying_saved_step':step,'complete_checkpoint_sha256':cp['sha256']})
                except Exception as exc:
                    failure=failure or {'reason':'FIELD_OR_CHECKPOINT_INVALID','detail':str(exc)};seen.add(step);request_stop('FAILURE',failure)
        if time.monotonic()-last_status>=15 or not alive:
            write_json(REPORT/'running.json',{'status':'RUNNING' if alive else 'EXITED','pid':proc.pid,
                'elapsed_s':time.monotonic()-start,'last_solve':{k:v for k,v in rows[-1].items() if k!='petsc_monitor'} if rows else None,
                'saved_steps':sorted(seen),'stop_request':stop_request,'failure':failure})
            last_status=time.monotonic()
        if not alive:break
        time.sleep(.5)
    code=proc.wait()
elapsed=time.monotonic()-start
history=parse_solver_log(log.read_text(errors='replace'),policy['dt_s'])
assert len(rows)==len(history['linear_solves'])
nl={r['step']:r for r in rows}
raw=resource.read_text();rss=re.search(r'Maximum resident set size \(kbytes\): (\d+)',raw)
clock=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)',raw)
wall=0
for part in clock[1].split(':'):wall=60*wall+float(part)
execution={'exit_code':code,'command':cmd,'PETSC_OPTIONS':opts['PETSC_OPTIONS'],'mpi_ranks':opts['mpi_ranks'],
    'OMP_NUM_THREADS':opts['OMP_NUM_THREADS'],'GNU_wall_s':wall,'elapsed_monotonic_s':elapsed,
    'peak_single_process_RSS_KiB':int(rss[1]),'memory_scope':'GNU time maximum individual process RSS, not MPI aggregate',
    'initial_step':initial_step,'last_step':max(nl),'steps_executed':len(nl),'linear':[{k:v for k,v in r.items() if k!='petsc_monitor'} for r in rows],
    'linear_statistics':stats([r['linear_iterations'] for r in rows]),'linear_failures':history['failed_linear_solves'],
    'nonlinear_failures':completed_nonlinear_failures(max(nl)),'ill_conditioned_warnings':history['ill_conditioned_warnings'],
    'stop_request':stop_request,'failure':failure,'log':str(log.relative_to(ROOT)),
    'resource_log':str(resource.relative_to(ROOT)),'full_monitors':'outputs/sv1_3/qc/cpu_linear_history.jsonl'}
write_json(REPORT/'cpu_execution.json',execution)
subprocess.run([sys.executable,'-B',str(ROOT/'scripts/sv13/validate_cpu.py')],check=True)

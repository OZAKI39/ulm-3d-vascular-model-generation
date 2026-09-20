"""One monitored, fresh native flow run; science verdict is produced from mirrored fields."""
import json, os, re, signal, subprocess, sys, tarfile, time
import xml.etree.ElementTree as ET
from pathlib import Path
from runner_remote import *
name=sys.argv[1];ranks=int(sys.argv[2]);backend=sys.argv[3];stack=sys.argv[4]
official=name.startswith('OFFICIAL_')
steady=name=='REAL_VASCULAR_GPU'
assert ranks in (1,4) and backend in ('CPU','GPU')
assert backend!='GPU' or ranks==1
assert re.fullmatch(r'[A-Za-z0-9_]+',name)
O=BASE/'outputs';O.mkdir(exist_ok=True)
manifest=json.loads((BASE/'configs/baseline_L_flow_input_manifest.json').read_text())
if not (O/'SV_MESH').exists():
    archive=BASE.parent/'sv1_3l/external/flow_inputs.tar.gz';assert digest(archive)==manifest['archive_sha256']
    with tarfile.open(archive) as t:t.extractall(O,filter='data')
    for f in manifest['files']:assert digest(O/f['destination'])==f['sha256']
case=O/name
if official:
    import shutil
    assert not case.exists();shutil.copytree(O/'official_fluid_gpu_smoke',case)
else:
    assert not case.exists();case.mkdir();(case/'solver.xml').write_bytes((O/'vascular_proof.xml').read_bytes())
if name in ('NEW_CPU_SHORT','NEW_GPU_SHORT'):
    short=json.loads((BASE/'configs/convergence_strategy.json').read_text())
    original=(case/'solver.xml').read_text()
    changed,n=re.subn(r'(<Number_of_time_steps>)\s*\d+\s*(</Number_of_time_steps>)',lambda m:m[1]+str(short['SHORT_REGRESSION_STEP'])+m[2],original)
    assert n==1
    (case/'solver.xml').write_text(changed)
if steady:
    assert backend=='GPU' and ranks==1
    strategy=json.loads((BASE/'configs/convergence_strategy.json').read_text());policy=strategy['steady_policy']
    original=(case/'solver.xml').read_text()
    for tag,value in [('Number_of_time_steps',policy['maximum_total_steps']),('Increment_in_saving_VTK_files',1)]:
        original,n=re.subn('(<'+tag+'>)\\s*\\d+\\s*(</'+tag+'>)',lambda m:m[1]+str(value)+m[2],original)
        assert n==1
    (case/'solver.xml').write_text(original)
tree=ET.parse(case/'solver.xml');dt=float(tree.findtext('.//Time_step_size'));steps=int(tree.findtext('.//Number_of_time_steps'))
tol=float(tree.findtext('.//Add_equation/Tolerance'));max_nl=int(tree.findtext('.//Add_equation/Max_iterations'))
base_options=json.loads((BASE/'configs/petsc_options.json').read_text())['PETSC_OPTIONS']
options=base_options+' -use_gpu_aware_mpi 0 -mat_type '+('aijcusparse' if backend=='GPU' else 'aij')+' -vec_type '+('cuda' if backend=='GPU' else 'standard')
profiling=name.startswith('PROFILE_')
benchmark=name.startswith('BENCH_')
if profiling:options+=' -log_view -log_view_gpu_time'
# Informational object views and true-residual monitors match the initial proof.
# Heavy profiling (log_view / GPU event timing / sampling) remains profile-only.
options+=' -vec_view ::ascii_info -mat_view ::ascii_info -ksp_view_mat ::ascii_info -ksp_view_rhs ::ascii_info -ksp_view_solution ::ascii_info'
if stack=='old':
    sv=json.loads((BASE/'configs/baseline_L_svmp_gpu_build.json').read_text());w=json.loads((BASE/'configs/baseline_L_compatibility_winner.json').read_text())
else:
    sv=load('svmp_'+stack+'_build');w=load('petsc_'+stack+'_build')
assert sv['status']=='PASS' and w['status']=='PASS'
mpi=json.loads((BASE/'configs/baseline_L_mpi_application_gate.json').read_text())
wrapper=w['candidate_wrapper']
rc=[Path.home()/'.petscrc',case/'.petscrc',case/'petscrc'];assert not any(p.exists() for p in rc),'AMBIENT_PETSC_OPTIONS'
runtime_env=env();runtime_env.update(LD_LIBRARY_PATH=sv['runtime_library_path'],PETSC_OPTIONS=options)
command=[wrapper,mpi['working_launcher'],'-n',str(ranks),wrapper,sv['executable'],'solver.xml']
if official and backend=='GPU':
    # Count actual PetscFinalize calls and inspect state at application MPI_Finalize.
    # This is one tiny official run, never a separate vascular profiling run.
    gdb=case/'lifecycle.gdb'
    gdb.write_text('set pagination off\nset confirm off\nset breakpoint pending on\nset $petsc_finalize_calls=0\nbreak PetscFinalize\ncommands\nsilent\nset $petsc_finalize_calls=$petsc_finalize_calls+1\ncontinue\nend\nbreak MPI_Finalize\nrun\nprintf "SV13N_PETSC_STATE initialized=%d finalized=%d count=%d\\n", *(unsigned char*)&PetscInitializeCalled, *(unsigned char*)&PetscFinalizeCalled, $petsc_finalize_calls\ncontinue\nquit $_exitcode\n')
    command=[wrapper,mpi['working_launcher'],'-n','1',wrapper,'gdb','--batch','-x',str(gdb),'--args',sv['executable'],'solver.xml']
log=L/(name+'.log');assert not log.exists()
from flow_parser import parse_solver_log
start_utc=time.time();start=time.monotonic();stop=None;offset=0;carry='';last_sample=0;gpu_samples=[]
with log.open('x') as f:
    process=subprocess.Popen(command,cwd=case,env=runtime_env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
    while process.poll() is None:
        time.sleep(.5)
        if steady and time.monotonic()-last_sample>=15:
            sample=subprocess.run(['nvidia-smi','--query-gpu=timestamp,memory.used,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=5)
            gpu_samples.append(dict(elapsed_s=time.monotonic()-start,exit_code=sample.returncode,raw=sample.stdout.strip()))
            write('vascular_gpu_device_snapshots',dict(samples=gpu_samples,scope='Device-wide inexpensive snapshots, not a profiling trace'))
            last_sample=time.monotonic()
        with log.open() as reader:reader.seek(offset);chunk=reader.read();offset=reader.tell()
        joined=carry+chunk
        if '\n' in joined:complete,carry=joined.rsplit('\n',1)
        else:complete='';carry=joined
        if stop is None:
            h=parse_solver_log(complete,dt)
            fatal=re.search(r'DIVERGED_\w+|PETSC ERROR|PC failed due to[^\n]*|(?<![a-z])(?:nan|inf)(?![a-z])',complete,re.I)
            nonlinear=[r for r in h['linear_solves'] if r['nonlinear_iteration']>=max_nl and (r['nonlinear_Ri_over_R0'] is None or r['nonlinear_Ri_over_R1'] is None or min(r['nonlinear_Ri_over_R0'],r['nonlinear_Ri_over_R1'])>tol)]
            heartbeat=case/'controller_heartbeat'
            monitor_lost=steady and time.monotonic()-start>120 and (not heartbeat.exists() or time.time()-heartbeat.stat().st_mtime>120)
            if fatal or h['failed_linear_solves'] or nonlinear or h['nonlinear_failure_messages'] or monitor_lost:
                stop={'reason':fatal.group() if fatal else 'STEADY_CONTROLLER_LOST' if monitor_lost else 'LINEAR_OR_NONLINEAR_FAILURE','wall_time_s':time.monotonic()-start}
                pending=case/'STOP_SIM.pending';pending.write_text('0\n');pending.replace(case/'STOP_SIM')
                if not steady:os.killpg(process.pid,signal.SIGTERM)
            elif not steady and time.monotonic()-start>14400:
                stop={'reason':'RUN_TIMEOUT_4H','wall_time_s':time.monotonic()-start};os.killpg(process.pid,signal.SIGTERM)
        elif not steady and time.monotonic()-start-stop['wall_time_s']>10:os.killpg(process.pid,signal.SIGKILL)
    code=process.wait()
elapsed=time.monotonic()-start;logtext=log.read_text(errors='replace');history=parse_solver_log(logtext,dt)
types=re.findall(r'(?:Mat|Vec|Matrix|Vector) Object:[^\n]*\n\s*type:\s*(\S+)',logtext)
mat_types=sorted(set(re.findall(r'(?:Mat|Matrix) Object:[^\n]*\n\s*type:\s*(\S+)',logtext)))
vec_types=sorted(set(re.findall(r'(?:Vec|Vector) Object:[^\n]*\n\s*type:\s*(\S+)',logtext)))
record={'status':'EXECUTED','name':name,'backend':backend,'mpi_ranks':ranks,'GPU_count':1 if backend=='GPU' else 0,'OMP_NUM_THREADS':1,
    'stack':stack,'solver_sha256':digest(sv['executable']),'PETSc_library_sha256':digest(Path(w['prefix'])/'lib/libpetsc.so'),'command':command,'cwd':str(case),'exit_code':code,'wall_time_s':elapsed,'log':'logs/'+log.name,'log_sha256':digest(log),
    'start_unix_s':start_utc,'end_unix_s':time.time(),'wall_clock_resolution_note':'Process completion detected at 0.5 s monitor interval; same policy for every candidate.',
    'PETSC_OPTIONS':options,'profiling':profiling,'benchmark':benchmark,'initial_state':'t=0','requested_steps':steps,
    'xml_sha256':digest(case/'solver.xml'),'monitor_stop':stop,'history':history,'matrix_types':mat_types,'vector_types':vec_types,
    'ambient_rc_files_checked_absent':[str(p) for p in rc],
    'results':[{'path':str(p.relative_to(BASE)),'sha256':digest(p),'size':p.stat().st_size} for p in sorted(case.glob('*-procs/result_*.vtu'))]}
write(name+'_execution',record)
print(json.dumps({'name':name,'exit_code':code,'wall_time_s':elapsed,'monitor_stop':stop,'mat_types':mat_types,'vec_types':vec_types,'result_count':len(record['results'])}),flush=True)
raise SystemExit(0 if code==0 and stop is None else 1)

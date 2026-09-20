"""One monitored, fresh native flow run; science verdict is produced from mirrored fields."""
import json, os, re, signal, subprocess, sys, tarfile, time
import xml.etree.ElementTree as ET
from pathlib import Path
from runner_remote import *
name=sys.argv[1];ranks=int(sys.argv[2]);backend=sys.argv[3]
assert ranks in (1,4) and backend in ('CPU','GPU')
assert backend!='GPU' or ranks==1
assert load('svmp_gpu_build')['status']=='PASS'
if name!='official_fluid_gpu_smoke':assert load('svmp_gpu_smoke')['status']=='PASS'
assert re.fullmatch(r'[A-Za-z0-9_]+',name)
O=BASE/'outputs';O.mkdir(exist_ok=True)
manifest=json.loads((BASE/'configs/flow_input_manifest.json').read_text())
if not (O/'SV_MESH').exists():
    archive=BASE/'external/flow_inputs.tar.gz';assert digest(archive)==manifest['archive_sha256']
    with tarfile.open(archive) as t:t.extractall(O,filter='data')
    for f in manifest['files']:assert digest(O/f['destination'])==f['sha256']
case=O/name
if name=='official_fluid_gpu_smoke':assert case.exists() and not list(case.glob('*-procs'))
else:
    assert not case.exists();case.mkdir();(case/'solver.xml').write_bytes((O/'vascular_proof.xml').read_bytes())
tree=ET.parse(case/'solver.xml');dt=float(tree.findtext('.//Time_step_size'));steps=int(tree.findtext('.//Number_of_time_steps'))
tol=float(tree.findtext('.//Add_equation/Tolerance'));max_nl=int(tree.findtext('.//Add_equation/Max_iterations'))
base_options=json.loads((BASE/'configs/petsc_options.json').read_text())['PETSC_OPTIONS']
options=base_options+' -use_gpu_aware_mpi 0 -mat_type '+('aijcusparse' if backend=='GPU' else 'aij')+' -vec_type '+('cuda' if backend=='GPU' else 'standard')
profiling=name.startswith('PROFILE_')
benchmark=name.startswith('BENCH_')
if profiling:options+=' -log_view -log_view_gpu_time'
if not benchmark:options+=' -vec_view ::ascii_info -mat_view ::ascii_info -ksp_view_mat ::ascii_info -ksp_view_rhs ::ascii_info -ksp_view_solution ::ascii_info'
if benchmark:
    for flag in ('-ksp_view','-options_view','-options_left'):options=options.replace(flag,'')
sv=load('svmp_gpu_build');w=load('compatibility_winner');mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text())
wrapper=w['candidate_wrapper']
rc=[Path.home()/'.petscrc',case/'.petscrc',case/'petscrc'];assert not any(p.exists() for p in rc),'AMBIENT_PETSC_OPTIONS'
runtime_env=env();runtime_env.update(LD_LIBRARY_PATH=sv['runtime_library_path'],PETSC_OPTIONS=options)
command=[wrapper,mpi['working_launcher'],'-n',str(ranks),wrapper,sv['executable'],'solver.xml']
log=L/(name+'.log');assert not log.exists()
from flow_parser import parse_solver_log
start=time.monotonic();stop=None;offset=0;carry=''
with log.open('x') as f:
    process=subprocess.Popen(command,cwd=case,env=runtime_env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
    while process.poll() is None:
        time.sleep(.5)
        with log.open() as reader:reader.seek(offset);chunk=reader.read();offset=reader.tell()
        joined=carry+chunk
        if '\n' in joined:complete,carry=joined.rsplit('\n',1)
        else:complete='';carry=joined
        if stop is None:
            h=parse_solver_log(complete,dt)
            fatal=re.search(r'DIVERGED_\w+|PETSC ERROR|PC failed due to[^\n]*|(?<![a-z])(?:nan|inf)(?![a-z])',complete,re.I)
            nonlinear=[r for r in h['linear_solves'] if r['nonlinear_iteration']>=max_nl and (r['nonlinear_Ri_over_R0'] is None or r['nonlinear_Ri_over_R1'] is None or min(r['nonlinear_Ri_over_R0'],r['nonlinear_Ri_over_R1'])>tol)]
            if fatal or h['failed_linear_solves'] or nonlinear or h['nonlinear_failure_messages']:
                stop={'reason':fatal.group() if fatal else 'LINEAR_OR_NONLINEAR_FAILURE','wall_time_s':time.monotonic()-start}
                (case/'STOP_SIM').write_text('0\n');os.killpg(process.pid,signal.SIGTERM)
            elif time.monotonic()-start>14400:
                stop={'reason':'RUN_TIMEOUT_4H','wall_time_s':time.monotonic()-start};os.killpg(process.pid,signal.SIGTERM)
        elif time.monotonic()-start-stop['wall_time_s']>10:os.killpg(process.pid,signal.SIGKILL)
    code=process.wait()
elapsed=time.monotonic()-start;logtext=log.read_text(errors='replace');history=parse_solver_log(logtext,dt)
types=re.findall(r'(?:Mat|Vec|Matrix|Vector) Object:[^\n]*\n\s*type:\s*(\S+)',logtext)
mat_types=sorted(set(re.findall(r'(?:Mat|Matrix) Object:[^\n]*\n\s*type:\s*(\S+)',logtext)))
vec_types=sorted(set(re.findall(r'(?:Vec|Vector) Object:[^\n]*\n\s*type:\s*(\S+)',logtext)))
record={'status':'EXECUTED','name':name,'backend':backend,'mpi_ranks':ranks,'GPU_count':1 if backend=='GPU' else 0,'OMP_NUM_THREADS':1,
    'command':command,'cwd':str(case),'exit_code':code,'wall_time_s':elapsed,'log':'logs/'+log.name,'log_sha256':digest(log),
    'PETSC_OPTIONS':options,'profiling':profiling,'benchmark':benchmark,'initial_state':'t=0','requested_steps':steps,
    'xml_sha256':digest(case/'solver.xml'),'monitor_stop':stop,'history':history,'matrix_types':mat_types,'vector_types':vec_types,
    'ambient_rc_files_checked_absent':[str(p) for p in rc],
    'results':[{'path':str(p.relative_to(BASE)),'sha256':digest(p),'size':p.stat().st_size} for p in sorted(case.glob('*-procs/result_*.vtu'))]}
write(name+'_execution',record)
print(json.dumps({'name':name,'exit_code':code,'wall_time_s':elapsed,'monitor_stop':stop,'mat_types':mat_types,'vec_types':vec_types,'result_count':len(record['results'])}),flush=True)
raise SystemExit(0 if code==0 and stop is None else 1)

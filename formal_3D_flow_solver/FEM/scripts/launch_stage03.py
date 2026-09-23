#!/usr/bin/env python3
"""Supervise the single MPI direct solve with live process-group RSS and cgroup limits."""
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import timestamp,write_json,command,sha256

base=ROOT/'outputs/stage03/reference'
for sub in ('solution','qc','metadata','checkpoints','logs'):(base/sub).mkdir(parents=True,exist_ok=True)
assert not (base/'metadata/production_attempt.json').exists(),'Only one formal factorization attempt is authorized'
pre=json.loads((ROOT/'outputs/stage03/preflight/assembly_verified.json').read_text());assert pre['status']=='PASS'

def memory():
    m={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if len(line.split())>=3}
    cg={}
    for name in ('memory.max','memory.current','memory.peak','memory.events','cpu.max'):
        path=Path('/sys/fs/cgroup')/name
        if path.exists():cg[name]=path.read_text().strip()
    available=m['MemAvailable']
    if cg.get('memory.max','max')!='max':available=min(available,int(cg['memory.max'])-int(cg['memory.current']))
    return {'timestamp':timestamp(),'meminfo_bytes':m,'cgroup':cg,'effective_available_bytes':available}

fresh=memory();required=pre['sum_rank_peak_rss_kib']*1024
initial=json.loads((ROOT/'outputs/stage03/resources/initial.json').read_text())
gate={'status':'PASS' if fresh['effective_available_bytes']>required else 'BLOCKED','fresh_memory':fresh,
      'assembled_process_peak_sum_bytes':required,'mumps_available':initial['mumps_available'],
      'hostname':platform.node(),'cpu':command(['lscpu']),'filesystem':command(['df','-B1',str(ROOT)]),
      'policy':'Available physical/cgroup memory must exceed measured assembly process memory; this is not a bound on LU fill-in. Monitor live process-group RSS and stop on exhausted cgroup headroom.',
      'gpu_used':False,'mpi_ranks':4,'omp_num_threads':1}
write_json(base/'metadata/resource_gate.json',gate)
assert gate['status']=='PASS' and gate['mumps_available']
cmd=[str(ROOT/'remote/.env/bin/mpiexec'),'-n','4',str(ROOT/'remote/.env/bin/python'),str(ROOT/'scripts/run_stage03.py'),'solve']
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
write_json(base/'metadata/production_attempt.json',{'timestamp':timestamp(),'command':cmd,'environment':{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},'attempt':1})
start=time.perf_counter();rss_peak=0;stopped=False;samples=[]
with (base/'logs/solve.stdout.txt').open('w') as out,(base/'logs/solve.stderr.txt').open('w') as err,(base/'logs/resources.jsonl').open('w') as log:
    proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=out,stderr=err,start_new_session=True)
    while proc.poll() is None:
        total=0;members=[]
        for p in Path('/proc').iterdir():
            if not p.name.isdigit():continue
            try:
                # /proc/PID/stat command names can contain spaces or parentheses.
                rest=(p/'stat').read_text().rsplit(')',1)[1].split()
                if int(rest[2])!=proc.pid:continue
                status=(p/'status').read_text().splitlines()
                rss=next((int(line.split()[1])*1024 for line in status if line.startswith('VmRSS:')),0)
                total+=rss;members.append({'pid':int(p.name),'rss_bytes':rss})
            except (OSError,ValueError,IndexError):continue
        rss_peak=max(rss_peak,total);now=memory()
        row={'elapsed_s':time.perf_counter()-start,'aggregate_process_group_rss_bytes':total,'members':members,
             'effective_available_bytes':now['effective_available_bytes'],'cgroup_memory_current':now['cgroup'].get('memory.current')}
        samples.append(row);log.write(json.dumps(row)+'\n');log.flush()
        if now['effective_available_bytes']<512*1024**2 and not stopped:
            stopped=True;os.killpg(proc.pid,signal.SIGTERM)
        time.sleep(.5)
    returncode=proc.wait()
text=(base/'logs/solve.stderr.txt').read_text()+(base/'logs/solve.stdout.txt').read_text()
resource_failed=stopped or returncode in (-9,137) or any(x in text.lower() for x in ('out of memory','cannot allocate memory','memory allocation','oom-kill','exit code: 9'))
result={'status':'PASS' if returncode==0 else ('BLOCKED' if resource_failed else 'FAIL'),
        'reason':'DIRECT_SOLVER_RESOURCE' if resource_failed else None,'returncode':returncode,
        'elapsed_time_s':time.perf_counter()-start,'peak_sampled_aggregate_rss_bytes':rss_peak,'sampling_interval_s':.5,
        'rusage_children_peak_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        'peak_definition':'Maximum sampled simultaneous RSS sum of the MPI process group; shared pages may be counted per process. Also retain child ru_maxrss (maximum individual child) and per-rank peaks.',
        'resource_guard_stopped':stopped,'fresh_memory_before':fresh,'memory_after':memory(),
        'stdout_sha256':sha256(base/'logs/solve.stdout.txt'),'stderr_sha256':sha256(base/'logs/solve.stderr.txt'),
        'gpu_used':False,'formal_attempts':1,'iterative_fallback':False}
write_json(base/'metadata/resources.json',result)
print(json.dumps(result,indent=2));print(text[-12000:])
sys.exit(returncode if returncode>0 else (1 if returncode<0 else 0))

"""Small adapter to the existing ledger/watchdog, with separate CPU accounting."""
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import threading
import time
from py_scripts.fluid_physics.common import PROJECT_ROOT, fingerprint, sha256_file, read_json, write_json, atomic_state, now, environment_identity
from py_scripts.fluid_physics.runner import run_attempt, shared_budget_state, exclusive_lock, group_members
from .physics import definition, lattice_definition


def paths(c):
    base=PROJECT_ROOT/'data/solver_benchmark'/c['campaign_id']
    runs=PROJECT_ROOT/'runs/solver_benchmark'/c['campaign_id']
    return base,runs


def plans(c):
    p,u=definition(c);env=environment_identity()
    dt=p['candidate']['dt_star']*u.t0
    gpu=dict(backend='Mirheo',method='SDPD',physics=p,dt_si=dt,
             steps=round(c['sampling']['formal_end_si']/dt),warmup_steps=round(c['sampling']['unforced_preparation_si']/dt),
             sample_steps=round(c['sampling']['interval_star']/p['candidate']['dt_star']),
             project_root=str(PROJECT_ROOT),source_commit=env['mirheo_commit'],binary_hash=env['library_sha256'],
             precision=env['precision'],worker_sha256=sha256_file(Path(__file__).with_name('mirheo_worker.py')),
             config_sha256=c['_config_sha256'])
    gpu['plan_sha256']=fingerprint(gpu)
    prior=paths(c)[0]/'frozen_plan.json'
    if prior.exists():
        approved=read_json(prior)['gpu_plan']
        skip={'worker_sha256','plan_sha256'}
        if {k:v for k,v in approved.items() if k not in skip}!={k:v for k,v in gpu.items() if k not in skip}:
            raise ValueError('SCIENTIFIC_PLAN_CHANGED')
        # Code corrections do not change the approved physics/budget. Record both hashes;
        # never replace the original plan, approval, or worker snapshot.
        gpu=approved
    tasks=[]
    for ranks in c['lbm']['ranks']:
        if ranks>len(os.sched_getaffinity(0)):continue
        for dx in c['lbm']['dx_candidates_si']:
            if ranks!=1 and dx!=c['lbm']['dx_candidates_si'][0]:continue
            l=lattice_definition(p,dx,c['lbm']['tau'])
            for rep in range(c['lbm']['repetitions']):
                task=dict(backend='HemoCell',method='D3Q19_Guo_BGK',role='main',rank_count=ranks,repetition=rep,
                          **l,steps=round(c['sampling']['formal_end_si']/l['dt_si']),
                          warmup_steps=round(c['sampling']['unforced_preparation_si']/l['dt_si']),
                          sample_steps=max(1,round(c['sampling']['interval_star']*u.t0/l['dt_si'])))
                task['task_id']=f'lbm_N{l["N"]}_r{ranks}_main_rep{rep+1}';tasks.append(task)
        l=lattice_definition(p,c['lbm']['dx_candidates_si'][0],c['lbm']['tau'])
        for rep in range(c['lbm']['repetitions']):
            tasks.append(dict(backend='HemoCell',method='D3Q19_Guo_BGK',role='low_io_cost',rank_count=ranks,repetition=rep,**l,
                steps=12000,warmup_steps=0,sample_steps=200,task_id=f'lbm_N{l["N"]}_r{ranks}_lowio_rep{rep+1}'))
    smoke={**tasks[0],'role':'smoke','steps':c['lbm']['smoke_steps'],'warmup_steps':0,'sample_steps':8,'task_id':'deployment_smoke'}
    return p,u,gpu,[smoke]+tasks


def gpu_budget(c,gpu):
    _,runs=paths(c);campaign=runs/'gpu'
    state=shared_budget_state(campaign,{'read_only_unregistered_scope':True,'shared_budget_pool':PROJECT_ROOT/c['shared_budget_pool']},task_id='sdpd_main')
    scope={'campaign_directory':str(campaign),'task_id':'sdpd_main','plan_sha256':gpu['plan_sha256']}
    matches=[r for r in state['authorization_records'] if r['scope']==scope]
    return dict(status='AUTHORIZED' if matches else 'PENDING_SCOPE_APPROVAL',authorization_scope=scope,
                requested_additional_seconds=c['budget']['gpu_requested_s'],single_task_limit_s=600,gpu_concurrency=1,
                shared_total_charged_or_reserved_s=state['total_charged_or_reserved_s'],
                original_remaining_s=state['original_remaining_s'],global_remaining_s=state['global_remaining_s'],
                new_benchmark_currently_authorized_s=sum(x['additional_seconds'] for x in matches),
                extension_accounts=state['extension_accounts'],
                reason='Old extensions are restricted to named restart/equilibration tasks. No permission is inferred from unallocated base balance; request one explicit 600 s new scope.',
                planned_steps=gpu['steps'],dt_si=gpu['dt_si'],planned_time_si=gpu['steps']*gpu['dt_si'],
                plan_sha256=gpu['plan_sha256'])


def freeze(c):
    base,runs=paths(c);base.mkdir(parents=True,exist_ok=True);runs.mkdir(parents=True,exist_ok=True)
    p,u,gpu,tasks=plans(c)
    frozen={'config_sha256':c['_config_sha256'],'physics':p,'gpu_plan':gpu,'cpu_tasks':tasks,'criteria':c['criteria'],'sampling':c['sampling']}
    path=base/'frozen_plan.json'
    if path.exists():
        if read_json(path)!=frozen:raise ValueError('FROZEN_PLAN_CHANGED; preserve prior campaign and inspect, never reset budget')
    else:write_json(path,frozen)
    request=gpu_budget(c,gpu)
    if not (base/'budget_request.json').exists():write_json(base/'budget_request.json',request)
    return frozen,request


def cpu_env():
    # Only this child receives system paths; parent and Mirheo environment untouched.
    return ['/usr/bin/env','-u','LD_LIBRARY_PATH','-u','PYTHONPATH','-u','PYTHONHOME','-u','HDF5_DIR',
            'PATH=/usr/bin:/bin','OMP_NUM_THREADS=1','OMPI_CC=/usr/bin/gcc','OMPI_CXX=/usr/bin/g++']


def xml_config(task,p):
    values={'N':task['N'],'bins':p['bins'],'dx':task['dx_si'],'dt':task['dt_si'],'nuP':p['nu_si'],'rhoP':p['rho_si'],'acceleration':p['acceleration_si']}
    domain=''.join(f'<{k}>{v!r}</{k}>' for k,v in values.items())
    return '<?xml version="1.0"?>\n<hemocell><parameters><outputDirectory>native_output</outputDirectory><logDirectory>log</logDirectory></parameters><domain>'+domain+'</domain><sim>'+f'<steps>{task["steps"]}</steps><warmupSteps>{task["warmup_steps"]}</warmupSteps><sampleEvery>{task["sample_steps"]}</sampleEvery></sim></hemocell>\n'


def process_tree_members(root_pid):
    # OpenMPI ranks can create distinct process groups: trace actual parentage.
    seen=set();queue=[root_pid]
    while queue:
        pid=queue.pop()
        if pid in seen:continue
        seen.add(pid)
        try:
            for thread in Path(f'/proc/{pid}/task').iterdir():
                queue.extend(int(x) for x in (thread/'children').read_text().split())
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    return sorted(seen)


def host_sample(pgid):
    members=process_tree_members(pgid);rows=[]
    for pid in members:
        try:
            status={line.split(':')[0]:line.split(':',1)[1].strip() for line in Path(f'/proc/{pid}/status').read_text().splitlines() if ':' in line}
            rank=None
            for part in Path(f'/proc/{pid}/environ').read_bytes().split(b'\0'):
                if part.startswith(b'OMPI_COMM_WORLD_RANK='):rank=int(part.split(b'=')[1])
            rows.append({'pid':pid,'name':status['Name'],'mpi_rank':rank,'rss_bytes':int(status.get('VmRSS','0 kB').split()[0])*1024,
                         'VmHWM_bytes':int(status.get('VmHWM','0 kB').split()[0])*1024,'cpus_allowed':status.get('Cpus_allowed_list')})
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    available=next(int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))
    return {'time_monotonic':time.perf_counter(),'processes':rows,'tree_rss_bytes':sum(x['rss_bytes'] for x in rows),'MemAvailable_bytes':available}


class HostMonitor:
    def __init__(self):self.stop=threading.Event();self.directory=None;self.samples=[]
    def run(self):
        while not self.stop.wait(.1):
            try:
                if self.directory and (self.directory/'process_identity.json').exists():
                    self.samples.append(host_sample(read_json(self.directory/'process_identity.json')['pgid']))
            except (OSError,ValueError):pass
    def __enter__(self):self.thread=threading.Thread(target=self.run,daemon=True);self.thread.start();return self
    def __exit__(self,*args):self.stop.set();self.thread.join(timeout=2)
    def record(self):
        ranks={}
        for s in self.samples:
            for p in s['processes']:
                if p['mpi_rank'] is not None:
                    key=str(p['mpi_rank']);prev=ranks.get(key,{'sampled_peak_rss_bytes':0})
                    ranks[key]={**p,'sampled_peak_rss_bytes':max(prev['sampled_peak_rss_bytes'],p['rss_bytes'])}
        return {'scope':'0.1 s sampled descendant process tree RSS sum across different process groups, includes shared pages in each process; not PSS or GPU VRAM',
                'coverage':'DESCENDANT_TREE',
                'sample_count':len(self.samples),'sampled_peak_tree_rss_bytes':max((x['tree_rss_bytes'] for x in self.samples),default=None),
                'minimum_WSL_available_bytes':min((x['MemAvailable_bytes'] for x in self.samples),default=None),
                'ranks':ranks,'samples':self.samples}


def runner_config(c,backend):
    return dict(campaign_id=c['campaign_id']+'_'+backend,
        budget=dict(campaign_limit_s=c['budget']['cpu_limit_s'] if backend=='cpu' else c['budget']['gpu_requested_s'],task_limit_s=600,
            minimum_launch_budget_s=5,graceful_margin_s=8,gpu_sample_interval_s=2,vram_reserve_MiB=1536,estimated_task_vram_MiB=512,host_reserve_MiB=1024),
        **({'shared_budget_pool':str(PROJECT_ROOT/c['shared_budget_pool']),'require_shared_budget':True} if backend=='gpu' else {}))


def seal(directory):
    hashes={str(p.relative_to(directory)):sha256_file(p) for p in directory.rglob('*') if p.is_file() and p.name!='output_sha256.json'}
    atomic_state(directory/'output_sha256.json',hashes)


def verify_run(directory):
    for name,sha in read_json(directory/'output_sha256.json').items():
        if sha256_file(directory/name)!=sha:raise ValueError('CACHE_OUTPUT_HASH_MISMATCH '+str(directory/name))


def execute_cpu(c,frozen):
    _,runs=paths(c);p=frozen['physics'];root=Path(c['hemocell_root']);binary=root/'build/benchmark/pure_fluid_benchmark'
    if not binary.is_file():raise FileNotFoundError('Run scripts/setup_hemocell.sh first')
    if subprocess.check_output(['git','-C',str(root/'vendor/HemoCell'),'rev-parse','HEAD'],text=True).strip()!=c['hemocell_commit']:raise ValueError('SOURCE_COMMIT_MISMATCH')
    results=[]
    for task in frozen['cpu_tasks']:
        ledger_path=runs/'cpu/budget_ledger.json'
        if ledger_path.exists() and sum(x.get('charged_s',x['reserved_s']) for x in read_json(ledger_path)['attempts'])>=c['budget']['cpu_limit_s']:break
        identity={'task':task,'config_sha256':c['_config_sha256'],'binary_sha256':sha256_file(binary)}
        with HostMonitor() as monitor:
            def command(d):
                monitor.directory=d
                write_json(d/'task.json',{**task,'binary_hash':sha256_file(binary),'source_commit':c['hemocell_commit']})
                (d/'config.xml').write_text(xml_config(task,p))
                cmd=cpu_env()+['/usr/bin/mpirun.openmpi','--bind-to','core','--map-by','core','--report-bindings','-np',str(task['rank_count']),str(binary),'config.xml']
                write_json(d/'command.json',cmd)
                return cmd
            # Same cross-solver lock also respects prior Mirheo campaigns.
            with exclusive_lock(PROJECT_ROOT/'runs/fluid_calibration/.gpu_exclusive.lock'):
                d,result,cached=run_attempt(runs/'cpu',runner_config(c,'cpu'),task['task_id'],identity,command,
                    c['budget']['cpu_task_allocation_s'],monitor_gpu=False)
        if d and not cached:write_json(d/'host_memory.json',monitor.record());seal(d)
        results.append({'task_id':task['task_id'],'directory':str(d) if d else None,'cached':cached,'execution':result})
        print(json.dumps({'task':task['task_id'],'status':result['status'],'cached':cached,'wall_s':result.get('elapsed_monotonic_s')},ensure_ascii=False),flush=True)
        if result['status']!='COMPLETED':break
    return results


def execute_gpu(c,frozen,*,explicit_retry=False):
    gpu=frozen['gpu_plan'];b=gpu_budget(c,gpu)
    if b['status']!='AUTHORIZED':return {'status':'NOT_RUN','reason':'PENDING_SCOPE_APPROVAL'}
    gpu={**gpu,'runtime_worker_sha256':sha256_file(Path(__file__).with_name('mirheo_worker.py')),
         'runtime_mpi_adapter_sha256':sha256_file(Path(__file__).with_name('native_mpi.py'))}
    _,runs=paths(c);campaign=runs/'gpu';pool_path=PROJECT_ROOT/c['shared_budget_pool']
    if explicit_retry:
        base=paths(c)[0];q=read_json(base/'explicit_retry_confirmation.json');request=read_json(base/'explicit_retry_request.json')
        if (q.get('authority')!='explicit_user_message' or q['user_message_sha256']!=sha256_file(Path(q['user_message_file']))
                or q['request_sha256']!=sha256_file(base/'explicit_retry_request.json') or q['maximum_retry_count']!=1
                or request['worker_sha256']!=gpu['runtime_worker_sha256'] or request['mpi_adapter_sha256']!=gpu['runtime_mpi_adapter_sha256']):
            raise ValueError('EXPLICIT_RETRY_EVIDENCE_MISMATCH')
        if len(read_json(campaign/'budget_ledger.json')['attempts'])!=1:
            raise ValueError('ONE_EXPLICIT_RETRY_ALREADY_USED')
    # Authorization records pre-exist and are validated by shared_budget_state.
    with exclusive_lock(PROJECT_ROOT/'runs/fluid_calibration/.gpu_exclusive.lock'):
        pool=read_json(pool_path)
        if not any(x['directory']==str(campaign) for x in pool['members']):
            pool['members'].append({'directory':str(campaign),'campaign_id':c['campaign_id']+'_gpu','ledger_required':False})
            atomic_state(pool_path,pool)
    with HostMonitor() as monitor:
        def command(d):
            monitor.directory=d;write_json(d/'task.json',gpu)
            worker=Path(__file__).with_name('mirheo_worker.py')
            (d/'mirheo_worker.py').write_bytes(worker.read_bytes())
            (d/'native_mpi.py').write_bytes(Path(__file__).with_name('native_mpi.py').read_bytes())
            launch='#!/bin/bash\nset -euo pipefail\nsource '+shlex.quote(str(PROJECT_ROOT/'scripts/activate_mirheo.sh'))+'\nexport SOLVER_BENCHMARK_AUTHORIZED='+shlex.quote(gpu['plan_sha256'])+'\nexec /usr/bin/mpirun.openmpi --bind-to none -np 2 '+shlex.quote(str(PROJECT_ROOT/'.venv/bin/python'))+' -B -u '+shlex.quote(str(d/'mirheo_worker.py'))+' --spec '+shlex.quote(str(d/'task.json'))+'\n'
            (d/'launch.sh').write_text(launch)
            return ['/bin/bash',str(d/'launch.sh')]
        # Failed startup is charged; explicit retry can use only this authorization's remainder.
        allocation=c['budget']['gpu_requested_s']
        if explicit_retry and (campaign/'budget_ledger.json').exists():
            allocation-=sum(x.get('charged_s',x['reserved_s']) for x in read_json(campaign/'budget_ledger.json')['attempts'])
        d,result,cached=run_attempt(campaign,runner_config(c,'gpu'),'sdpd_main',gpu,command,allocation,strict_allocation=True,retry_failed=explicit_retry)
    if d and not cached:write_json(d/'host_memory.json',monitor.record());seal(d)
    return {'directory':str(d) if d else None,'cached':cached,**result}


def register_authorization(c,confirmation_path):
    """Administrative entrypoint after explicit real user authorization only."""
    frozen,_=freeze(c);base,runs=paths(c);confirmation_path=Path(confirmation_path).resolve();q=read_json(confirmation_path)
    request=read_json(base/'budget_request.json');evidence=Path(q['user_message_file']).resolve();message=evidence.read_text().strip()
    if q.get('authority')!='explicit_user_message' or not message or q['additional_seconds']!=c['budget']['gpu_requested_s']:
        raise ValueError('EXPLICIT_USER_MESSAGE_AND_EXACT_SECONDS_REQUIRED')
    if q['budget_request_sha256']!=sha256_file(base/'budget_request.json') or q['scope']!=request['authorization_scope']:
        raise ValueError('AUTHORIZATION_SCOPE_OR_PLAN_MISMATCH')
    record={'schema_version':1,'id':fingerprint(q),'authority':'explicit_user_message','user_message_file':str(evidence),
            'user_message_sha256':sha256_file(evidence),'user_message':message,'additional_seconds':q['additional_seconds'],
            'scope':q['scope'],'single_task_limit_s':600,'gpu_concurrency':1,'budget_request_file':str(base/'budget_request.json'),
            'budget_request_sha256':q['budget_request_sha256'],'confirmation_file':str(confirmation_path),
            'confirmation_sha256':sha256_file(confirmation_path),'recorded_at':now()}
    pool_path=PROJECT_ROOT/c['shared_budget_pool']
    with exclusive_lock(PROJECT_ROOT/'runs/fluid_calibration/.gpu_exclusive.lock'):
        pool=read_json(pool_path)
        from py_scripts.fluid_physics.budget_authorizations import authorization_records
        if any(x['id']==record['id'] for x in authorization_records(pool)):raise ValueError('DUPLICATE_AUTHORIZATION')
        dest=runs/'gpu/authorizations'/record['id'];dest.mkdir(parents=True,exist_ok=False)
        write_json(dest/'shared_pool_before.json',pool);write_json(dest/'authorization.json',record)
        pool.setdefault('authorization_records',[]).append({'path':str(dest/'authorization.json'),'sha256':sha256_file(dest/'authorization.json')})
        atomic_state(pool_path,pool)
    return str(dest/'authorization.json')

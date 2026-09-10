"""Small adapter around the existing watchdog; no implicit authorization."""
import json
from pathlib import Path
import shutil
import time
import yaml
from py_scripts.fluid_physics.common import sha256_file, fingerprint, write_json, atomic_state, now
from py_scripts.fluid_physics.runner import run_attempt
from py_scripts.single_rbc_benchmark.workflow import child_env
from .forensics import ROOT, read

DEFAULT_CONFIG = ROOT/'py_scripts/single_rbc_benchmark_repaired.yaml'


def load_config(path=DEFAULT_CONFIG):
    path = Path(path).resolve()
    c = yaml.safe_load(path.read_text())
    name = c['campaign_id']
    if not name or name in ('.','..') or '/' in name or '\\' in name:
        raise ValueError('INVALID_CAMPAIGN')
    if c['repair']['legacy_campaign_id'] == name:
        raise ValueError('LEGACY_CAMPAIGN_IS_READ_ONLY')
    if any(c['budget'][k] for k in ('gpu_authorized_s','cpu_solver_authorized_s','native_build_authorized_s')):
        raise ValueError('CONFIG_IS_NOT_AN_AUTHORIZATION')
    c['_config_path'], c['_config_sha256'] = str(path), sha256_file(path)
    return c


def paths(c):
    return tuple(ROOT / folder / 'single_rbc_repair' / c['campaign_id'] for folder in ('data','runs','test_code/outputs'))


def verify_protection(c):
    b,_,_ = paths(c)
    manifest = read(b/'protection_before.json')['files']
    changed = [p for p,h in manifest.items() if not Path(p).is_file() or sha256_file(p)!=h]
    result = dict(status='PASS' if not changed else 'FAILED', checked_files=len(manifest), changed=changed, checked_at=now())
    atomic_state(b/'protection_current.json',result)
    if changed:
        raise ValueError('PROTECTED_FILES_CHANGED')
    return result


def require_authorization(c):
    b,_,_ = paths(c)
    p = b/'authorization.json'
    if not p.exists():
        raise RuntimeError('NEW_REPAIR_AUTHORIZATION_REQUIRED: old 788.296746 seconds are scoped to the archived configuration')
    a = read(p)
    request = read(b/'authorization_request.json')
    if a.get('approved') is not True or a.get('request_sha256')!=fingerprint(request) or a.get('scope')!='single_rbc_repair_ABC':
        raise RuntimeError('REPAIR_AUTHORIZATION_SCOPE_MISMATCH')
    if not a.get('user_quote') or not a.get('recorded_at'):
        raise RuntimeError('USER_APPROVAL_EVIDENCE_REQUIRED')
    plan = read(b/'frozen_benchmark_plan.json')
    for p,h in plan['executable_artifact_sha256'].items():
        if not Path(p).is_file() or sha256_file(p)!=h:
            raise RuntimeError('FROZEN_EXECUTABLE_CHANGED: '+p)
    if c['_config_sha256']!=plan['config_sha256']:
        raise RuntimeError('FROZEN_CONFIG_CHANGED')
    return a


def counted_run(c, task, category, command_builder, allocation, identity):
    require_authorization(c)
    b,runs,_ = paths(c)
    limits = {'gpu':c['budget']['gpu_proposal_s'],'build':c['budget']['native_build_proposal_s'], 'cpu':c['budget']['cpu_solver_proposal_s']}
    budget = dict(campaign_limit_s=limits[category], task_limit_s=c['budget']['max_solver_task_s'],
                  minimum_launch_budget_s=3, graceful_margin_s=10, gpu_sample_interval_s=3,
                  vram_reserve_MiB=1536, estimated_task_vram_MiB=1800, host_reserve_MiB=1024)
    runner = dict(campaign_id=c['campaign_id']+'_'+category,budget=budget)
    return run_attempt(runs/category,runner,task,identity,command_builder,allocation,
                       monitor_gpu=category=='gpu',retry_failed=False,strict_allocation=True)


def isolated_build(c):
    require_authorization(c)
    b,_,_ = paths(c)
    record = read(b/'isolated_build.json')
    if record.get('status')=='BUILT':
        p = Path(record['library'])
        if not p.is_file() or sha256_file(p)!=record['library_sha256']:
            raise RuntimeError('ISOLATED_LIBRARY_CHANGED')
        return p
    src,build = b/'native/source',b/'native/build'
    args = ['/usr/bin/cmake','-S',str(src),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release',
            '-DCMAKE_C_COMPILER=/usr/bin/gcc-12','-DCMAKE_CXX_COMPILER=/usr/bin/g++-12',
            '-DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.6/bin/nvcc','-DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-12',
            '-DCMAKE_CUDA_ARCHITECTURES=89','-DMIR_CUDA_ARCH_NAME=8.9','-DMIR_ENABLE_STACKTRACE=OFF',
            '-DMIR_BUILD_TESTS=OFF','-DMIR_DOUBLE_PRECISION=OFF','-DMIR_MEMBRANE_DOUBLE=OFF',
            '-DPython_EXECUTABLE='+str(ROOT/'.venv/bin/python')]
    for task,cmd,seconds in [('configure_isolated',args,120),('compile_isolated',['/usr/bin/cmake','--build',str(build),'--parallel','2'],1680)]:
        def create(d,cmd=cmd):
            write_json(d/'command.json',cmd)
            return child_env()+cmd
        _,execution,_ = counted_run(c,task,'build',create,seconds,dict(command=cmd,patch_sha256=sha256_file(b/'native/local_before_halo_v2.patch')))
        if execution['status']!='COMPLETED':
            raise RuntimeError('ISOLATED_BUILD_FAILED_NO_AUTOMATIC_RETRY: '+task)
    candidates = list(build.rglob('libmirheo.cpython-312-*.so'))
    if len(candidates)!=1:
        raise RuntimeError('AMBIGUOUS_BUILT_LIBRARY')
    library = candidates[0]
    atomic_state(b/'isolated_build.json',{**record,'status':'BUILT','built_at':now(),'library':str(library),
                                         'library_sha256':sha256_file(library),'installed_library_replaced':False})
    return library


def short_spec(c, task, library=None):
    b,_,_ = paths(c)
    spec = dict(config=c,role='main',dt=.001,seed=20260910,mesh=str(b/'common_reference.off'),
                prep_steps=5000,steps=0,sample_steps=100,stop_after_relaxation=True,
                bouncer_policy=task['bouncer_policy'],debug_level=task.get('debug_level',1),diagnostic_state=True,
                ks=c['mirheo_membrane']['ks_candidate'],kb=c['mirheo_membrane']['kb_candidate'])
    if library:
        spec.update(native_library=str(library),native_library_sha256=sha256_file(library))
    spec['plan_sha256'] = fingerprint(read(b/'frozen_benchmark_plan.json'))
    return spec


def short_diagnostics(c):
    require_authorization(c)
    verify_protection(c)
    b,_,_ = paths(c)
    plan = read(b/'frozen_benchmark_plan.json')
    results = []
    for task in plan['first_executable_queue']:
        library = isolated_build(c) if task['library']=='isolated' else None
        spec = short_spec(c,task,library)
        identity = dict(spec=spec,worker_sha256=sha256_file(Path(__file__).parent/'mirheo_worker.py'))
        def create(d):
            start = time.perf_counter()
            write_json(d/'spec.json',spec)
            write_json(d/'case_creation.json',{'elapsed_s':time.perf_counter()-start,'new_inputs':True})
            shutil.copyfile(Path(__file__).parent/'mirheo_worker.py',d/'worker_source.py')
            return child_env(gpu=True)+['PYTHONPATH='+str(ROOT),'SINGLE_RBC_AUTHORIZED='+spec['plan_sha256'],
                    'RBC_REPAIR_NATIVE_TRACE='+('1' if task.get('trace') else '0'),
                    '/usr/bin/mpirun.openmpi','--bind-to','none','-np','2',str(ROOT/'.venv/bin/python'),'-B',
                    '-m','py_scripts.single_rbc_repair.mirheo_worker','--spec',str(d/'spec.json')]
        directory,execution,cached = counted_run(c,task['id'],'gpu',create,task['allocation_s'],identity)
        results.append(dict(task=task['id'],directory=str(directory),execution=execution,status=execution['status'],
                            cached=cached,new_repeat=not cached))
        atomic_state(b/'short_diagnostics_execution.json',{'runs':results,'qualification':'none; CPU stored-state review required before any longer task'})
        if execution['status']!='COMPLETED':
            raise RuntimeError('SHORT_DIAGNOSTIC_STOPPED; distinct tests require review of the preserved failure, no automatic retry')
    return dict(status='SHORT_RUNS_RECORDED_REQUIRE_GATE_A_REVIEW',runs=results,
                warning='A short result never establishes a full-range fix. Material matching and Gamma=4 runs remain gated.')


def qualification(runs, material, screen, target=4.):
    if material!='PASS' or screen!='PASS':
        return None
    for solver in ('HemoCell','Mirheo'):
        selected = [r for r in runs if r['solver']==solver and r.get('completed') and not r.get('cached')
                    and r.get('strain_end')==target]
        if len({r['run_id'] for r in selected}) < 2:
            return None
    return 'ELIGIBLE_FOR_MEASURED_TIMING_RATIO'

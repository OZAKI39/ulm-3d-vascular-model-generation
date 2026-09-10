"""Sequential native tasks. Requires the new campaign's explicit authorization."""
from pathlib import Path
import copy
import json
import shutil
from py_scripts.fluid_physics.common import PROJECT_ROOT, read_json, write_json, sha256_file, fingerprint
from .physics import paths, sample_steps
from .workflow import run, run_cpu, child_env
from .analysis import viscosity, affine_response

def gpu_task(c,name,role,allocation,*,strict=False,parameters=None,runtime_env=None):
    base,_,_=paths(c);dt=c['dpd']['strict_dt' if strict else 'dt']
    spec=dict(config=c,role=role,dt=dt,seed=20260910+(1 if name.endswith('_2') else 0),mesh=str(base/'common_reference.off'),
              **sample_steps(c,dt))
    spec.update(parameters or {})
    if runtime_env:
        if runtime_env!={'CUDA_DEVICE_MAX_CONNECTIONS':'1','CUDA_DEVICE_MAX_COPY_CONNECTIONS':'1'}:raise ValueError('UNREVIEWED_RUNTIME_DIAGNOSTIC')
        spec['runtime_environment']=runtime_env
    spec['worker_sha256']=sha256_file(Path(__file__).with_name('mirheo_worker.py'));spec['plan_sha256']=fingerprint(spec)
    def create(d):
        write_json(d/'spec.json',spec)
        snapshot=d/'source';snapshot.mkdir()
        for filename in ('mirheo_worker.py','physics.py','MIRHEO_LICENSE'):
            shutil.copyfile(Path(__file__).with_name(filename),snapshot/filename)
        return child_env(gpu=True)+[k+'='+v for k,v in (runtime_env or {}).items()]+['SINGLE_RBC_AUTHORIZED='+spec['plan_sha256'],'PYTHONPATH='+str(PROJECT_ROOT),
            '/usr/bin/mpirun.openmpi','--bind-to','core','--map-by','core','-np','2',str(PROJECT_ROOT/'.venv/bin/python'),'-B','-m','py_scripts.single_rbc_benchmark.mirheo_worker','--spec',str(d/'spec.json')]
    return run(c,name,spec,create,allocation,gpu=True)

def execute(c):
    """Run the frozen sequential list; verified cached evidence is never a repeat."""
    from py_scripts.solver_benchmark.workflow import verify_run
    base,runs,_=paths(c)
    auth=read_json(base/'authorization.json')
    if auth.get('approved') is not True or auth.get('scope')!='single_rbc_DPD_shear' or auth.get('config_sha256')!=c['_config_sha256']:
        raise RuntimeError('AUTHORIZATION_SCOPE_MISMATCH')
    for name in ('common_reference.off','screen_frozen.json','preflight.json','formal_frozen.json','fluid_measured.json','membrane_measured.json'):
        if not (base/name).exists():raise RuntimeError('FROZEN_PREPARATION_REQUIRED '+name)
    frozen=read_json(base/'formal_frozen.json')
    if frozen['config_sha256']!=c['_config_sha256']:raise ValueError('FORMAL_PLAN_CHANGED')
    nu=frozen['nu'];params=frozen['membrane'];tasks=[
        ('cpu',dict(id='empty_hemocell',cell=False,ranks=2,allocation_s=80)),
        ('gpu',dict(id='empty_dpd',role='empty',allocation_s=550)),
        ('cpu',dict(id='main_hemocell_1',ranks=2,allocation_s=120)),
        ('gpu',dict(id='main_mirheo_1',role='main',allocation_s=950)),
        ('cpu',dict(id='main_hemocell_2',ranks=2,allocation_s=120)),
        ('gpu',dict(id='main_mirheo_2',role='main',allocation_s=950)),
        ('cpu',dict(id='strict_hemocell',ranks=2,strict=True,allocation_s=220)),
        ('gpu',dict(id='strict_mirheo',role='main',strict=True,allocation_s=950))]
    results=[]
    for backend,t in tasks:
        directory=runs/'solver'/t['id']
        if (directory/'execution.json').exists():
            verify_run(directory);r=read_json(directory/'execution.json')
            results.append(dict(directory=str(directory),status=r['status'],cached=True))
            print(t['id'],r['status'],'existing sealed evidence, not a new cold repeat',flush=True)
            continue
        ledger=read_json(runs/'solver/budget_ledger.json')
        charged=sum(a.get('charged_s',a.get('reserved_s',0)) for a in ledger['attempts'])
        remaining=c['budget']['proposed_total_solver_s']-charged
        if remaining<10:
            results.append(dict(task=t['id'],status='NOT_RUN_BUDGET_LIMIT'));break
        t={**t,'allocation_s':min(t['allocation_s'],remaining)}
        print('START',t['id'],'allocation_s',t['allocation_s'],flush=True)
        if backend=='cpu':d,r,cached=run_cpu(c,t,nu)
        else:
            parameters=dict(params)
            if t['id']=='strict_mirheo':
                parameters['steps']=round(frozen['strict_mirheo_endpoint_strain']/c['protocol']['shear_rate']/c['dpd']['strict_dt'])
            d,r,cached=gpu_task(c,t['id'],t['role'],t['allocation_s'],strict=t.get('strict',False),parameters=parameters)
        results.append(dict(directory=str(d),status=r['status'],cached=cached))
        print('FINISH',t['id'],r['status'],r['elapsed_monotonic_s'],flush=True)
        # No automatic retry. Each failure remains charged, and independent tasks continue.
    if not (base/'execution_batch.json').exists():write_json(base/'execution_batch.json',dict(tasks=results,automatic_retries=0))
    from .reporting import export
    return export(c)

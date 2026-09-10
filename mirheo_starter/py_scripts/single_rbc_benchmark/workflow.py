"""Small campaign adapter; reuse the existing watchdog, ledger and RSS monitor."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from py_scripts.fluid_physics.common import PROJECT_ROOT, fingerprint, sha256_file, read_json, write_json, atomic_state, environment_identity, resource_snapshot, now
from py_scripts.fluid_physics.runner import run_attempt
from py_scripts.solver_benchmark.workflow import HostMonitor, seal, verify_run
from .physics import paths, lattice, xml_case, rbc_files, units, geometry, read_off, groups

def environment(c):
    h=Path(c['hemocell']['root']);m=environment_identity()
    git={}
    for root in (h/'vendor/HemoCell',PROJECT_ROOT/'vendor/Mirheo'):
        git[str(root)]={}
        for key,args in [('commit',['rev-parse','HEAD']),('status',['status','--short'])]:
            r=subprocess.run(['git','-C',str(root),*args],capture_output=True,text=True)
            git[str(root)][key]=dict(returncode=r.returncode,stdout=r.stdout.strip(),stderr=r.stderr.strip())
    source=read_json(h/'metadata/sources.json');archive=h/'downloads'/('palabos-'+source['palabos_commit']+'.tar.gz')
    archive_sha=sha256_file(archive)
    if archive_sha!=source['archive_sha256']:raise ValueError('PALABOS_SOURCE_ARCHIVE_HASH_MISMATCH')
    palabos=dict(kind='archive plus the existing HemoCell patch; not a separate Git checkout',commit=source['palabos_commit'],archive_sha256=archive_sha,manifest=str(h/'metadata/sources.json'),patch_sha256=sha256_file(h/'vendor/HemoCell/patch/palabos.patch'))
    return dict(recorded_at=now(),resources=resource_snapshot(),cpu_affinity=sorted(os.sched_getaffinity(0)),mirheo=m,git=git,palabos=palabos,
                hemocell_library_sha256=sha256_file(h/'build/upstream/libhemocell.a'),threads=dict(OMP=1,OPENBLAS=1,MKL=1),
                precision={'HemoCell':'double, verify native completion sizeof(T)','Mirheo':m['precision']})

def child_env(gpu=False):
    env=['/usr/bin/env','-u','LD_LIBRARY_PATH','-u','PYTHONPATH','-u','PYTHONHOME','-u','HDF5_DIR',
         'PATH=/usr/bin:/bin','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1','PYTHONNOUSERSITE=1']
    if gpu:env+=['LD_LIBRARY_PATH=/usr/local/cuda-12.6/lib64:/usr/lib/wsl/lib']
    return env

def run(c,task_id,identity,make_command,allocation,*,gpu=False,preparation=False,retry=False):
    base,runs,_=paths(c)
    if gpu:
        auth=base/'authorization.json'
        if not auth.exists():raise RuntimeError('NEW_GPU_AUTHORIZATION_REQUIRED; old budgets are read only')
        a=read_json(auth)
        if a.get('approved') is not True or a.get('config_sha256')!=c['_config_sha256'] or a.get('scope')!='single_rbc_DPD_shear':raise RuntimeError('AUTHORIZATION_SCOPE_MISMATCH')
    campaign=runs/('preparation' if preparation else 'solver')
    budget=dict(campaign_limit_s=360 if preparation else c['budget']['proposed_total_solver_s'],task_limit_s=c['budget']['proposed_max_task_s'],
                minimum_launch_budget_s=3,graceful_margin_s=8,gpu_sample_interval_s=3,vram_reserve_MiB=1536,estimated_task_vram_MiB=1024,host_reserve_MiB=1024)
    config=dict(campaign_id=c['campaign_id']+('_preparation' if preparation else '_solver'),budget=budget)
    # This ledger is deliberately not enrolled into any older shared budget pool.
    identity={**identity,'backend':'GPU' if gpu else 'CPU','config_sha256':c['_config_sha256']}
    creation_s=0.
    with HostMonitor() as monitor:
        def create(d):
            nonlocal creation_s
            monitor.directory=d;t=time.perf_counter()
            try:return make_command(d)
            finally:creation_s=time.perf_counter()-t
        d,record,cached=run_attempt(campaign,config,task_id,identity,create,allocation,monitor_gpu=gpu,retry_failed=retry,strict_allocation=True)
    if d and not cached:
        write_json(d/'case_creation.json',dict(wall_s=creation_s,scope='fresh input/configuration files and command construction before solver process; added to E2E, outside solver budget clock'))
        write_json(d/'host_resources.json',monitor.record());seal(d)
    return d,record,cached

def compile_case(c):
    root=Path(c['hemocell']['root']);src=root/'cases/single_rbc_shear_benchmark';build=root/'build/single_rbc_shear_benchmark'
    identity={p.name:sha256_file(p) for p in src.iterdir() if p.is_file()}
    def configure(d):return child_env()+['/usr/bin/cmake','-S',str(src),'-B',str(build),'-DCMAKE_C_COMPILER=/usr/bin/gcc','-DCMAKE_CXX_COMPILER=/usr/bin/g++']
    key=fingerprint(identity)[:8]
    d,r,_=run(c,'cmake_configure_'+key,identity,configure,30,preparation=True)
    if r['status']!='COMPLETED':raise RuntimeError('CMAKE_CONFIGURE_FAILED '+str(d))
    d,r,_=run(c,'compile_case_'+key,identity,lambda d:child_env()+['/usr/bin/cmake','--build',str(build),'--parallel','1'],c['budget']['one_time_cpu_compile_limit_s'],preparation=True)
    if r['status']!='COMPLETED':raise RuntimeError('CASE_COMPILE_FAILED '+str(d))
    return build/'single_rbc_shear_benchmark'

def run_cpu(c,task,nu=None):
    binary=Path(c['hemocell']['root'])/'build/single_rbc_shear_benchmark/single_rbc_shear_benchmark'
    nu=nu if nu is not None else c['material_check']['nu_candidate']
    xml=xml_case(c,nu,strict=task.get('strict',False),cell=task.get('cell',True),probe=task.get('probe',False),short_steps=task.get('short_steps'))
    identity=dict(task=task,binary_sha256=sha256_file(binary),xml_sha256=fingerprint(xml),material=c['hemocell']['material'],resolved_config=c)
    def create(d):
        (d/'config.xml').write_text(xml);rbc_files(c,d);write_json(d/'spec.json',identity)
        snapshot=d/'source';snapshot.mkdir()
        for source in (Path(c['hemocell']['root'])/'cases/single_rbc_shear_benchmark').iterdir():
            if source.is_file():shutil.copyfile(source,snapshot/source.name)
        return child_env()+['/usr/bin/mpirun.openmpi','--bind-to','core','--map-by','core','-np',str(task.get('ranks',2)),str(binary),'config.xml']
    return run(c,task['id'],identity,create,task['allocation_s'])

def protection(c):
    base,_,_=paths(c);before=read_json(base/'protection_before.json');changed=[]
    for path,sha in before['files'].items():
        p=Path(path)
        if not p.is_file() or sha256_file(p)!=sha:changed.append(path)
    record=dict(recorded_at=now(),status='PASS' if not changed else 'FAILED',files_checked=len(before['files']),changed=changed)
    # New campaign audit only; no old record is modified.
    atomic_state(base/'protection_current.json',record);return record

def preflight(c,prepare_cpu=False):
    base,runs,review=paths(c)
    for d in (base,runs,review):d.mkdir(parents=True,exist_ok=True)
    if not (base/'environment.json').exists():write_json(base/'environment.json',environment(c))
    if prepare_cpu:
        binary=Path(c['hemocell']['root'])/'build/single_rbc_shear_benchmark/single_rbc_shear_benchmark'
        compile_case(c)
        for task in [dict(id='candidate1_native_response',ranks=1,probe=True,short_steps=0,allocation_s=20),
                     dict(id='candidate1_r2_1024',ranks=2,short_steps=1024,allocation_s=35),
                     dict(id='candidate1_r1_256',ranks=1,short_steps=256,allocation_s=15)]:
            d,r,_=run_cpu(c,task)
            if r['status']!='COMPLETED':raise RuntimeError('CPU_PREFLIGHT_FAILED '+str(d))
    mesh_run=runs/'solver/candidate1_native_response/reference.off';g=None
    if mesh_run.exists():
        if not (base/'common_reference.off').exists():shutil.copyfile(mesh_run,base/'common_reference.off')
        g=geometry(*read_off(base/'common_reference.off'))
        if not (base/'geometry.json').exists():write_json(base/'geometry.json',g)
    outcomes=[]
    for p in runs.rglob('execution.json'):
        r=read_json(p);outcomes.append(dict(directory=str(p.parent),elapsed_s=r['elapsed_monotonic_s'],status=r['status'],completion=read_json(p.parent/'completion.json') if (p.parent/'completion.json').exists() else None))
    if not (base/'screen_frozen.json').exists():write_json(base/'screen_frozen.json',dict(config_sha256=c['_config_sha256'],criteria=c['screen'],protocol=c['protocol'],freeze_stage='before any GPU or formal timing',nu_definition='measured DPD material result admissible [0.4,4.0], fixed before formal runs'))
    record=dict(stage='CPU_PREFLIGHT_GPU_NOT_AUTHORIZED',config_sha256=c['_config_sha256'],geometry=g,units=units(),
                fluid=dict(status='MEASUREMENT_PENDING',nu_candidate_only=c['material_check']['nu_candidate'],method='native DPD double periodic Poiseuille force profile, independently repeated at half dt'),
                membrane=dict(status='PARTIAL',reason='affine native force and independent relaxation checks required; HO bending reference differs from constant-angle Kantor; membrane inertia differs'),
                proposed_budget=c['budget'],protocol=c['protocol'],screen=c['screen'],tasks=outcomes,
                qualified_speedup=None,gpu_authorized=False,protection=protection(c),
                estimates=dict(cpu_main_each_s=[75,100],cpu_strict_s=[145,160],gpu_main_each_s=[300,550],gpu_strict_s=[600,900],
                    gpu_estimate_status='UNMEASURED: based on particle/step count and previous native runs; no claim of verified throughput',
                    common_fluid_particles_estimated=8*24**3,membrane_vertices=642,lattice_nodes=48*48*49,
                    priorities=['independent material and coupling checks','both full cold workflows','second cold repeat','own stricter dt','independent analysis'],
                    budget_definition='sum of sequential solver process wall seconds; CPU seconds are not multiplied by rank count; compile separately measured; failed preflights charged'))
    if not (base/'preflight.json').exists():write_json(base/'preflight.json',record)
    if (base/'authorization.json').exists():
        a=read_json(base/'authorization.json');record['gpu_authorized']=a.get('approved') is True and a.get('config_sha256')==c['_config_sha256']
        record['stage']='AUTHORIZED_CAMPAIGN_OBSERVATION'
    if (base/'fluid_measured.json').exists():
        record['fluid']={'status':'MEASURED','evidence':read_json(base/'fluid_measured.json')}
    if (base/'membrane_measured.json').exists():record['membrane']=read_json(base/'membrane_measured.json')
    if (base/'native_failure_final.json').exists():
        record['stage']='COMPUTATION_STOPPED_FINAL_EVIDENCE'
        record['final_outcome']=read_json(base/'native_failure_final.json')
    atomic_state(base/'preflight_latest.json',record);return record

"""Full-RBC provenance and preflight. Never imports Mirheo or starts a solver."""
from common import *
import math, shutil, subprocess,struct

EXPECTED_LIBRARY = 'd45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee'
BUILD_ID = 'mirheo-sm120-20731713865ae510'
PURPOSE = 'One cold full RBC verification: preparation then Gamma=4; one GPU; no retry'
TOOL_NAMES = ['common.py','remote_rbc.py','rbc_full_worker.py','rbc_checks.py','rbc_h5.py',
              'rbc_geometry.py','rbc_geometry_sources.json','cloud_geometry.py','rbc_report.py']


def validate_protocol(s):
    c=s['config'];dt=s['dt'];p=c['protocol']
    if not (s.get('continuous_observation') and s['role']=='main' and s['bouncer_policy']=='shared'):
        raise ValueError('OLD_OR_INCOMPATIBLE_WORKER_PROTOCOL')
    if (dt,p['shear_rate'],p['end_strain'],p['relaxation_time']) != (.0005,.02,4.,30.):
        raise ValueError('LATEST_REPAIR_PARAMETERS_CHANGED')
    shear=round(p['end_strain']/p['shear_rate']/dt);prep=round(p['relaxation_time']/dt)
    assert shear==400000 and prep==60000 and c['dpd']['wall_relax_steps']==2000
    assert s['ks']==3 and s['kb']==8 and c['mirheo_membrane']['belonging_correction_every']==0
    native_dt=struct.unpack('f',struct.pack('f',dt))[0]
    return dict(dt=dt,native_float32_dt=native_dt,native_shear_time_star=shear*native_dt,native_strain_target_step=shear*native_dt*p['shear_rate'],prep_steps=prep,shear_steps=shear,wall_hidden_steps=4000,
                prep_time_star=prep*dt,shear_time_star=shear*dt,shear_rate=p['shear_rate'],target_strain=4,
                preparation_sample_steps=1000,shear_sample_steps=10000,maximum_membrane_frames=102,
                full_fluid_frames='only the two native phase-initial states; no full-fluid trajectory',
                membership_scope='up to 192 sampled IDs at phase initial and returned endpoints; not strict impermeability')


REMOTE_PROBE = '''
from pathlib import Path
import json,sys,hashlib,os,shutil,subprocess,socket
r=json.loads(sys.argv[1]);b=Path(r['build_dir']);i=read(b/'evidence/compile_identity.json')
assert read(b/'build_status.json')['status']=='CLOUD_BUILD_PASS'
assert i['runtime_binary_hash']==r['expected_library'] and sha(i['runtime_binary_path'])==r['expected_library'],'RUNTIME_BINARY_MISMATCH'
assert i['original_source_sha256']==r['original_source_sha256'] and i['target_architectures']==['120'] and i['precision_bits']==32,'BUILD_SOURCE_MISMATCH'
for entry in r['source_files']:
 assert sha(Path(r['resolved_cloud_source'])/entry['relative'])==entry['sha256'],'UPLOADED_SOURCE_CHANGED'
for entry in r['python_inputs']:
 assert sha(entry['remote'])==entry['sha256'],'UPLOADED_PYTHON_INPUT_CHANGED'
for rel,h in r['protected_physics_files'].items():
 assert sha(Path(r['work_source'])/rel)==h,'BUILD_PHYSICS_SOURCE_MISMATCH'
def snapshot(root):
 rows=[]
 for parent,dirs,names in os.walk(root,followlinks=False):
  for name in sorted(dirs+names):
   p=Path(parent)/name;s=p.lstat();row=dict(path=p.relative_to(root).as_posix(),mode=s.st_mode,mtime_ns=s.st_mtime_ns)
   if p.is_symlink():row['link']=os.readlink(p)
   elif p.is_file():row.update(size=s.st_size,sha256=sha(p))
   rows.append(row)
 return dict(files=len(rows),fingerprint=object_hash(sorted(rows,key=lambda x:x['path'])))
gpu=subprocess.run(['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total,memory.used,memory.free','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=15)
apps=subprocess.run(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],capture_output=True,text=True,timeout=15)
assert gpu.returncode==apps.returncode==0
own=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  argv=(p/'cmdline').read_bytes().replace(b'\\0',b' ').decode(errors='replace')
  if '/workspace/bloodflow/cloud_runs/' in argv and any(x in argv for x in ['worker.py','remote_smoke.py','remote_rbc.py']):own.append(dict(pid=int(p.name),command=argv))
 except (FileNotFoundError,PermissionError,ProcessLookupError):pass
cgroup={name:(Path('/sys/fs/cgroup')/name).read_text().strip() for name in ['memory.max','memory.current','cpu.max','cpuset.cpus.effective']}
commands={}
for name,argv in [('python',[r['python'],'-B','--version']),('mpi',['/usr/bin/mpirun','--version'])]:
 q=subprocess.run(argv,capture_output=True,text=True,timeout=10);commands[name]=dict(argv=argv,exit_code=q.returncode,stdout=q.stdout)
print(json.dumps(dict(checked_at=now(),hostname=socket.gethostname(),uid=os.getuid(),identity=i,gpu=gpu.stdout.strip(),gpu_compute_apps=apps.stdout.strip(),disk=shutil.disk_usage('/workspace')._asdict(),cgroup=cgroup,own_active_jobs=own,commands=commands,project_venv=snapshot(Path('/workspace/bloodflow/.venv')),smoke_budget=read('/workspace/bloodflow/metadata/cloud_smoke_budget_20260910.json'))))
'''


def preflight():
    import cloud_run as cr
    cr.local_environment();base=cr.BASE;directory=base/'rbc_full';directory.mkdir(exist_ok=True)
    build=read(base/'build_plan.json');manifest=read(cr.UPLOAD/'transfer_manifest.json')
    native=read(cr.UPLOAD/'metadata/native_source_locations.json');mapping=read(cr.UPLOAD/'metadata/path_mapping.json')
    tree,remote,key,project=cr.resolve_source(native,mapping,manifest)
    assert build['build_id']==BUILD_ID and build['original_source_path']==tree['path']
    rows=[x for x in manifest['files'] if x['path'].startswith(key+'/')]
    changed=[];source_files=[]
    for x in rows:
        rel=x['path'][len(key)+1:];p=Path(tree['path'])/rel
        if not p.is_file() or sha(p)!=x['sha256'] or sha(cr.UPLOAD/'source'/x['path'])!=x['sha256']: changed.append(rel)
        source_files.append(dict(relative=rel,sha256=x['sha256']))
    expected={x['relative'] for x in source_files}
    for parent,dirs,names in os.walk(tree['path']):
        dirs[:]=[x for x in dirs if x!='.git']
        for name in names:
            if name=='.git': continue
            rel=(Path(parent)/name).relative_to(tree['path']).as_posix()
            if rel not in expected:
                # This zero-byte placeholder was explicitly excluded as build-aux in the frozen transfer.
                excluded='units/extern/googletest/googlemock/build-aux/.keep'
                if rel==excluded and (Path(parent)/name).stat().st_size==0:
                    evidence=read(cr.UPLOAD/'metadata/excluded_files.json')['exclusions']
                    assert any(e['source_path']==str((Path(parent)/name).parent) and e['reason']=='Default build directory exclusion: build-aux' for e in evidence)
                else:changed.append('NEW:'+rel)
    if changed: raise ValueError('BUILD_SOURCE_MISMATCH:'+repr(changed))
    digest=object_hash([{k:x[k] for k in ['path','size','sha256','execute_bits']} for x in rows])
    assert digest==build['original_source_sha256']
    project_map=[x for x in mapping['roots'] if tree['path'].startswith(x['local']+'/')][0]
    original_spec=Path(native['latest_spec']['path']);s=read(original_spec)
    receipt=Path(native['explicit_delivery_receipt']['path']);receipt_data=read(receipt)
    files_to_check=[original_spec,receipt,Path(s['mesh'])]
    local_project=Path(project_map['local'])
    for suffix in ['py_scripts/single_rbc_repair/continuous_worker.py','py_scripts/single_rbc_repair/continuous_protocol.py',
                   'py_scripts/single_rbc_repair/coupling.py','py_scripts/solver_benchmark/native_mpi.py',
                   'py_scripts/single_rbc_repair/geometry_checks.py','py_scripts/single_rbc_repair/preparation_quality.py',
                   'py_scripts/single_rbc_benchmark/physics.py','py_scripts/single_rbc_benchmark/quality.py',
                   'data/single_rbc_repair/rbc_repair_20260910T131105Z/material_matching.json']:
        files_to_check.append(local_project/suffix)
    inputs=[]
    for p in files_to_check:
        rel=p.relative_to(local_project);frozen=cr.UPLOAD/project_map['staging']/rel
        if sha(p)!=sha(frozen): raise ValueError('PYTHON_OR_INPUT_UPDATE_REQUIRES_EXPLICIT_WHITELIST:'+str(p))
        inputs.append(dict(local=str(p),remote=project_map['remote']+'/'+rel.as_posix(),sha256=sha(p)))
    material=read(files_to_check[-1]);counts=validate_protocol(s)
    worker=base/'rbc_full_worker.py'
    assert 'CLOUD_FULL_RBC_V1' in worker.read_text() and 'continuous_protocol import evolve_stage' in worker.read_text(),'OLD_WORKER_REJECTED'
    tools={name:sha(base/name) for name in TOOL_NAMES}
    spec=json.loads(json.dumps(s));spec.pop('plan_sha256',None)
    spec.update(mesh=project_map['remote']+'/'+Path(s['mesh']).relative_to(local_project).as_posix(),steps=counts['shear_steps'],sample_steps=counts['shear_sample_steps'],
                prep_sample_steps=1000,stop_after_relaxation=False,cloud_full_protocol='CLOUD_FULL_RBC_V1',preparation_criteria=material['common_preparation_criteria'],
                native_library='/workspace/bloodflow/.venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so',native_library_sha256=EXPECTED_LIBRARY)
    spec['config']['protocol']['cold_repetitions']=1
    # Compute the common area radius without importing project code or numerical libraries.
    lines=Path(s['mesh']).read_text().splitlines();nv,nf,_=map(int,lines[1].split());v=[list(map(float,row.split())) for row in lines[2:2+nv]];area=0
    for row in lines[2+nv:2+nv+nf]:
        _,a,b,c=map(int,row.split());e=[v[b][k]-v[a][k] for k in range(3)];f=[v[c][k]-v[a][k] for k in range(3)]
        cross=[e[1]*f[2]-e[2]*f[1],e[2]*f[0]-e[0]*f[2],e[0]*f[1]-e[1]*f[0]];area+=math.sqrt(sum(q*q for q in cross))/2
    spec['reference_radius']=math.sqrt(area/(4*math.pi))
    smoke=read('/home/lzy/projects/cloud_results/cloud-smoke-20260910T221827Z/C_repair/completion.json')
    naive=smoke['relaxation_s']/250*(counts['prep_steps']+counts['shear_steps'])+smoke['setup_s']
    plan=dict(protocol='CLOUD_FULL_RBC_V1',build_id=BUILD_ID,original_source_sha256=digest,original_source_path=tree['path'],resolved_cloud_source=remote,
              spec=spec,python=build['python'],cloud_project_root=project,tool_hashes=tools,python_inputs=inputs,
              controller_hashes={name:sha(base/name) for name in ['cloud_run.py','rbc_plan.py','prepare_full_worker.py']},
              source_files=source_files,counts=counts,expected_output_bytes=512*1024**2,
              estimate=dict(smoke_extrapolated_process_s=naive,planning_range_s=[360,1500],basis='250 RBC steps of archived cloud smoke only; added CPU frame screening/native trace; unvalidated long-run estimate; not a completion guarantee'),
              limits=dict(process_seconds=1800,concurrency=1,formal_attempts=1,min_preflight_free_bytes=5*1024**3,min_running_free_bytes=3*1024**3),
              authority=dict(status='PENDING_SEPARATE_FULL_RUN_APPROVAL',purpose=PURPOSE,prior='cloud_smoke_budget_20260910 only authorizes small tests; local RBC budgets excluded'),
              criteria=dict(preparation=material['common_preparation_criteria'],numerical=material['criteria'],
                            hard_stop=['native CUDA/MPI/collision error','nonfinite state','proper self intersection','membrane wall crossing','degenerate triangle area<=1e-12','WLC extension>=1','fluid population change','time/disk/memory protection'],
                            soft_policy='A/V >2% records FAILED screen but preparation continues to its fixed cap; enter shear only if preparation windows and geometry pass. During shear minor screen exceedance is recorded and does not alter parameters.',
                            incomplete_preparation='stop at t*=30; no shear; no extra preparation or retry'),
              handoff=dict(strategy='existing two sequential coordinators in the same cold process',preserved=['all fluid positions/velocities sorted by old ID','wall positions; prescribed velocities changed to +/-0.24','membrane coordinates/velocities','original reference mesh'],
                           regenerated=['particle IDs via FromArray','initial member classification','membrane oldPositions=current initial coordinates','forces/task graph/cell lists','native stochastic generator state'],
                           seamless_state_continuity=False,checkpoint_loaded=False),
              scientific_baseline=dict(status=receipt_data.get('status'),physical_model_validation='NOT_MATCHED; convergence NOT_TESTED',qualified_speedup=None),
              environment=dict(HWLOC_COMPONENTS='-opencl,-gl',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',RBC_REPAIR_NATIVE_TRACE='1'),
              environment_evidence='Archived smoke and real MPI checks required -opencl,-gl; -opencl alone hung in GL discovery; no new CUDA concurrency restriction')
    plan['plan_sha256']=object_hash(plan)
    request=dict(build,expected_library=EXPECTED_LIBRARY,source_files=source_files,python_inputs=inputs)
    probe=cr.ssh_code((base/'common.py').read_text()+REMOTE_PROBE,request)
    blocks=[]
    if probe['own_active_jobs']: blocks.append('OWN_ACTIVE_JOB')
    if probe['gpu_compute_apps']: blocks.append('GPU_COMPUTE_ALREADY_ACTIVE')
    if probe['disk']['free']-plan['expected_output_bytes']<5*1024**3: blocks.append('DISK_HEADROOM')
    if '5090' not in probe['gpu']: blocks.append('GPU_IDENTITY_CHANGED')
    prior_gpu=read(base/'records/preflight.json')['commands']['gpu']['stdout'].split(',')[1].strip()
    if probe['gpu'].split(',')[1].strip()!=prior_gpu:blocks.append('GPU_UUID_CHANGED')
    if probe['uid']!=0:blocks.append('REMOTE_USER_CHANGED')
    if int(probe['gpu'].split(',')[-1])<4096: blocks.append('GPU_MEMORY_HEADROOM')
    cg=probe['cgroup']
    if cg['memory.max']!='max' and int(cg['memory.max'])-int(cg['memory.current'])<2*1024**3: blocks.append('CGROUP_MEMORY_HEADROOM')
    report=dict(status='PREFLIGHT_PASS_EXECUTION_PENDING_AUTHORIZATION' if not blocks else 'PREFLIGHT_BLOCKED',checked_at=now(),blockers=blocks,
                plan_sha256=plan['plan_sha256'],native_files_verified=len(rows),python_inputs_verified=len(inputs),source_differences=[],cloud=probe,
                prior_smoke_evidence=validate_result('/home/lzy/projects/cloud_results/cloud-smoke-20260910T221827Z'))
    version=directory/'plans'/(plan['plan_sha256']+'.json')
    if not version.exists():write(version,plan)
    write(directory/'frozen_plan.json',plan);write(directory/'preflight.json',report)
    print(json.dumps(dict(status=report['status'],plan=str(directory/'frozen_plan.json'),plan_sha256=plan['plan_sha256'],counts=counts,estimate=plan['estimate'],output_bytes=plan['expected_output_bytes'],cloud_free_bytes=probe['disk']['free'],blockers=blocks),ensure_ascii=False),flush=True)
    return plan,report


def authorization(base, plan):
    p=Path(base)/'rbc_full/authorization.json'
    if not p.is_file(): raise ValueError('FULL_RBC_AUTHORIZATION_REQUIRED: one new 1800-second full cloud task; smoke/local balances excluded')
    a=read(p)
    if a.get('plan_sha256')!=plan['plan_sha256'] or a.get('purpose')!=PURPOSE or a.get('process_seconds')!=1800 or a.get('attempts')!=1 or not a.get('user_approval_quote'):
        raise ValueError('FULL_RBC_AUTHORIZATION_SCOPE_MISMATCH')
    if (Path(base)/'rbc_full/launch_receipt.json').exists(): raise ValueError('FULL_RBC_ATTEMPT_ALREADY_RESERVED_USE_STATUS_FETCH')
    return a

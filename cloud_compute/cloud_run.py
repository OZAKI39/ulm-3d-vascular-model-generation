#!/usr/bin/env python3
"""Local SSH/rsync control of source-bound builds, smoke and one authorized full RBC task."""
from common import *
import argparse,subprocess,shlex,sys,shutil,re,difflib,platform,getpass
BASE=Path('/home/lzy/projects/cloud_compute')
UPLOAD=Path('/home/lzy/projects/cloud_upload/20260910T212439Z')
SSH=['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=15','-o','StrictHostKeyChecking=yes','-o','ClearAllForwardings=yes','-o','RemoteCommand=none','vast-mirheo']

def ssh_code(code,args=None,timeout=240):
    prefix='import sys\nsys.argv=["remote",'+repr(json.dumps(args or {}))+']\n'
    r=subprocess.run(SSH+[shlex.join(['python3','-B','-'])],input=(prefix+code).encode(),capture_output=True,timeout=timeout)
    if r.returncode:raise RuntimeError('SSH_REMOTE_FAILURE '+r.stderr.decode(errors='replace')[-1500:]+' '+r.stdout.decode(errors='replace')[-2000:])
    return json.loads(r.stdout)
def local_environment():
    assert getpass.getuser()=='lzy' and os.environ.get('HOME')=='/home/lzy' and 'wsl' in platform.release().lower(),'LOCAL_WSL_REQUIRED'
    r=subprocess.run(['ssh','-G','vast-mirheo'],capture_output=True,text=True,check=True)
    actual={}
    for line in r.stdout.splitlines():
        a=line.split(None,1)
        if len(a)==2 and a[0] in ['hostname','port','user','identityfile']:actual.setdefault(a[0],[]).append(a[1])
    assert actual['hostname']==['204.111.105.196'] and actual['port']==['50271'] and actual['user']==['root'],'SSH_TARGET_CHANGED'
    assert actual['identityfile'] in [['~/.ssh/id_ed25519_vast'],['/home/lzy/.ssh/id_ed25519_vast']],'SSH_IDENTITY_CHANGED'
    return actual
def resolve_source(native,mapping,manifest):
    candidates=[x for x in native['native_trees'] if x['path']==native['latest_source']]
    if len(candidates)!=1 or native['status']!='CONFIRMED_LATEST_NATIVE_SOURCE':raise ValueError('SOURCE_AMBIGUITY')
    item=candidates[0]
    if '/vendor/Mirheo' in item['path'] or '/native/source' not in item['path']:raise ValueError('OLD_VENDOR_SOURCE_REJECTED')
    maps=[x for x in mapping['roots'] if item['path'].startswith(x['local']+'/')]
    if len(maps)!=1:raise ValueError('PATH_MAPPING_MISSING')
    m=maps[0];relative=Path(item['path']).relative_to(m['local']).as_posix()
    key=Path(m['staging']).relative_to('source').as_posix()+'/'+relative
    if not any(x['path'].startswith(key+'/') for x in manifest['files']):raise ValueError('SOURCE_NOT_IN_MANIFEST')
    return item,m['remote']+'/'+relative,key,m['remote']

def preflight(revise=False):
    local=local_environment();probe=ssh_code((BASE/'remote_preflight.py').read_text());write(BASE/'records/preflight.json',probe)
    manifest=read(UPLOAD/'transfer_manifest.json');native=read(UPLOAD/'metadata/native_source_locations.json');mapping=read(UPLOAD/'metadata/path_mapping.json')
    receipt=read('/home/lzy/projects/cloud_upload/runtime_logs/20260910T212439Z-final-delivery.json')
    assert receipt['sha256sum_exit_code']==0 and receipt['status']=='PASS'
    tree,resolved,key,cloud_project=resolve_source(native,mapping,manifest)
    source=UPLOAD/'source'/key;rows=[x for x in manifest['files'] if x['path'].startswith(key+'/')]
    for x in rows:assert sha(UPLOAD/'source'/x['path'])==x['sha256'],'LOCAL_FROZEN_SOURCE_CHANGED'
    digest=object_hash([{k:x[k] for k in ['path','size','sha256','execute_bits']} for x in rows])
    original_version=(source/'mirheo/version.py').read_text();version=re.search(r'mir_version\s*=\s*"([^"]+)"',original_version)[1]
    edits=[];patches=[]
    for relative in ['cmake/version.cmake','setup.py','src/CMakeLists.txt']:
        old=(source/relative).read_text()
        if relative=='cmake/version.cmake':
            new='''# Cloud archive adapter: explicit provenance, no fabricated Git repository.
macro(getMirheoVersion version_str version_number)
  if(NOT DEFINED MIR_ARCHIVED_VERSION)
    message(FATAL_ERROR "MIR_ARCHIVED_VERSION is required for this Git-free archive")
  endif()
  set(${version_str} "${MIR_ARCHIVED_VERSION}")
  string(REGEX REPLACE "^v" "" ${version_number} "${MIR_ARCHIVED_VERSION}")
endmacro()
macro(getMirheoSHA1 sha1)
  if(NOT DEFINED MIR_ARCHIVED_REVISION)
    message(FATAL_ERROR "MIR_ARCHIVED_REVISION is required for archive provenance")
  endif()
  set(${sha1} "${MIR_ARCHIVED_REVISION}")
endmacro()
'''
        elif relative=='setup.py':
            needle="glob.glob(ext.sourcedir + '/build/src/mirheo/bindings/libmirheo.cpython*.so')"
            assert old.count(needle)==1
            new=old.replace(needle,"glob.glob(os.environ.get('MIRHEO_PREBUILT_LIBRARY', ext.sourcedir + '/build/src/mirheo/bindings/libmirheo.cpython*.so'))")
        else:
            needle='cuda_select_nvcc_arch_flags(BUGGED_ARCH_FLAGS ${MIR_CUDA_ARCH_NAME})'
            assert old.count(needle)==1
            replacement='''# CMake 3.28 FindCUDA accepts only single-digit architecture major versions.
# Keep its original path for old GPUs; explicitly guard and target Blackwell 12.0.
if(MIR_CUDA_ARCH_NAME STREQUAL "12.0")
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.8)
    message(FATAL_ERROR "sm_120 requires CUDA Toolkit 12.8 or newer")
  endif()
  set(BUGGED_ARCH_FLAGS "-gencode;arch=compute_120,code=sm_120")
else()
  cuda_select_nvcc_arch_flags(BUGGED_ARCH_FLAGS ${MIR_CUDA_ARCH_NAME})
endif()'''
            new=old.replace(needle,replacement)
        edits.append(dict(path=relative,original_sha256=sha(source/relative),replacement=new,replacement_sha256=hashlib.sha256(new.encode()).hexdigest()))
        patches.append(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+relative,tofile='b/'+relative)))
    options=['-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CUDA_ARCHITECTURES=120','-DMIR_CUDA_ARCH_NAME=12.0','-DMIR_BUILD_TESTS=OFF',
        '-DMIR_ENABLE_STACKTRACE=OFF','-DMIR_DOUBLE_PRECISION=OFF','-DMIR_MEMBRANE_DOUBLE=OFF','-DMIR_ROD_DOUBLE=OFF',
        '-DMIR_ENABLE_LTO=OFF','-DMIR_USE_NVTX=OFF','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON','-DPYBIND11_FINDPYTHON=ON',
        '-DPython_EXECUTABLE=/workspace/bloodflow/.venv/bin/python','-DPython_INCLUDE_DIR=/usr/include/python3.12',
        '-DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc','-DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++',
        '-DCMAKE_CXX_COMPILER=/usr/bin/g++','-DCMAKE_C_COMPILER=/usr/bin/gcc','-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda',
        '-DMPI_C_COMPILER=/usr/bin/mpicc','-DMPI_CXX_COMPILER=/usr/bin/mpicxx','-DHDF5_C_COMPILER_EXECUTABLE=/usr/bin/h5pcc',
        '-DHDF5_PREFER_PARALLEL=ON','-DMIR_ARCHIVED_VERSION='+version,'-DMIR_ARCHIVED_REVISION='+tree['head']+'-dirty-cloud']
    toolchain={k:probe['commands'][k]['stdout'] for k in ['nvcc','gcc','gxx','cmake','mpi','python','hdf5_link']}
    fingerprint=dict(source_sha256=digest,toolchain=toolchain,options=options,compatibility_patch=''.join(patches),gpu_architecture='sm_120',dependencies=['numpy==1.26.4'])
    build_id='mirheo-sm120-'+object_hash(fingerprint)[:16];build='/workspace/bloodflow/build/'+build_id;work='/workspace/bloodflow/work/'+build_id+'/source'
    protection=[]
    for parent,dirs,names in os.walk(UPLOAD,followlinks=False):
        for name in dirs+names:
            p=Path(parent)/name;relative=p.relative_to(UPLOAD).as_posix()
            if p.is_symlink():protection.append(dict(path=relative,type='symlink',target=os.readlink(p)))
            elif p.is_file():protection.append(dict(path=relative,sha256=sha(p),execute_bits=p.stat().st_mode&0o111))
    physics={x['path'][len(key)+1:]:x['sha256'] for x in rows if x['path'][len(key)+1:].startswith('src/') and not x['path'].endswith(('version.cpp','CMakeLists.txt'))}
    plan=dict(build_id=build_id,created_at=now(),original_source_path=tree['path'],resolved_cloud_source=resolved,original_source_sha256=digest,
        source_revision=tree['head'],source_manifest_sha256=sha(UPLOAD/'transfer_manifest.json'),upload_dir='/workspace/bloodflow/uploads/20260910T212439Z',
        cloud_project_root=cloud_project,python='/workspace/bloodflow/.venv/bin/python',build_dir=build,work_source=work,
        configure_argv=['cmake','-S',work,'-B',build,'-G','Ninja',*options],compatibility_edits=edits,compatibility_patch=''.join(patches),
        protected_physics_files=physics,upload_protection=sorted(protection,key=lambda x:x['path']),python_dependencies=['numpy==1.26.4'],
        toolchain=toolchain,build_fingerprint=object_hash(fingerprint),estimated_build_growth_bytes=3*1024**3,
        protected_tool_binaries=['/usr/local/cuda/bin/nvcc','/usr/bin/gcc','/usr/bin/g++','/usr/bin/cmake','/usr/bin/mpirun'],
        limits=dict(build_seconds=1800,test_seconds=60,total_test_seconds=180,concurrency=1,build_parallelism=2,min_run_free_bytes=3*1024**3),
        authority='Current pasted user deployment instruction; new cloud operation guard, no historical RBC balance reused',
        smoke_plan=dict(A='two-rank hostname and real Allreduce/Bcast/Barrier; installed native import and CUDA synchronization',B=dict(kind='periodic DPD liquid',domain=[4,4,4],steps=100,parameters_from='latest A6 spec.config.dpd'),C=dict(kind='short continuous deformable RBC preparation',optional=True,steps=250,dt=.0005,wall_steps_per_wall=2000,parameters_from='latest A6 spec; same membrane, bounce-back, mass and liquid parameters',scope='runtime path only; no scientific calibration or Gamma=4')),
        baseline_science=dict(status='GATE_A_UNRESOLVED_HALF_DT_PREPARATION_COMPLETED_QUALITY_FAILED',material_match='NOT_MATCHED',same_quality_same_endpoint_complete=False))
    assert probe['uid']==0 and probe['python_environment']['exists'] and not probe['python_environment']['symlink']
    assert 'sm_120' in probe['commands']['architectures']['stdout'] and 'release 12.8' in toolchain['nvcc']
    assert 'Parallel HDF5: yes' in probe['commands']['hdf5']['stdout']
    assert not probe['commands']['gpu_processes']['stdout'].strip(),'GPU already has a compute task; do not launch'
    assert disk_allowed(probe['disk']['free'],'build',plan['estimated_build_growth_bytes']),'INSUFFICIENT_BUILD_SPACE'
    path=BASE/'build_plan.json'
    if path.exists() and read(path)['build_id']!=build_id:
        assert revise,'Different existing build plan; inspect before replacing'
        prior=read(path);record=ssh_code('from pathlib import Path\nimport json,sys\nr=json.loads(sys.argv[1]);print((Path(r["build"])/"build_status.json").read_text())',dict(build=prior['build_dir']))
        assert record['status']=='FAILED','Do not replace a running/successful build plan'
        write(BASE/'records'/('build_plan_'+prior['build_id']+'.json'),prior)
        plan['supersedes_failed_configuration']=prior['build_id'];plan['work_source_reuse_from']=prior['work_source']
        plan['prior_deployment_elapsed_s']=record.get('cumulative_elapsed_s',record['elapsed_s'])
        write(path,plan)
    elif not path.exists():write(path,plan)
    print(json.dumps(dict(status='PREFLIGHT_PASS',build_id=build_id,source=resolved,source_sha256=digest,gpu=probe['commands']['gpu']['stdout'].strip(),disk=probe['disk'],limits=plan['limits'],smoke_plan=plan['smoke_plan'],plan=str(path)),ensure_ascii=False),flush=True)

def deploy(resume=False):
    local_environment();plan=read(BASE/'build_plan.json');build=plan['build_dir'];identity=dict(build_id=plan['build_id'],build_fingerprint=plan['build_fingerprint'])
    response=ssh_code('''from pathlib import Path
import json,sys,os,hashlib
r=json.loads(sys.argv[1]);p=Path(r['build']);assert p.parent==Path('/workspace/bloodflow/build') and not p.is_symlink()
identity=p/'identity.json'
if p.exists():
 assert identity.is_file() and json.loads(identity.read_text())==r['identity'],'BUILD_ID_COLLISION'
 status=json.loads((p/'build_status.json').read_text()) if (p/'build_status.json').is_file() else {'status':'PREPARED'}
 if status['status']=='CLOUD_BUILD_PASS':
  i=json.loads((Path(status['evidence_directory'])/'compile_identity.json').read_text())
  assert hashlib.sha256(Path(i['runtime_binary_path']).read_bytes()).hexdigest()==i['runtime_binary_hash'],'INSTALLED_BUILD_CHANGED'
  status['reused_after_binary_hash_verification']=True
 print(json.dumps(status))
else:
 p.mkdir(parents=True);(p/'tools').mkdir();(p/'evidence').mkdir();identity.write_text(json.dumps(r['identity']))
 print(json.dumps({'status':'NEW'}))
''',dict(build=build,identity=identity))
    attempt=1
    if response['status']=='FAILED' and resume:
        attempt=response.get('attempt',1)+1
        assert attempt<=3 and response.get('cumulative_elapsed_s',response.get('elapsed_s',0))<plan['limits']['build_seconds'],'BUILD_RETRY_OR_TIME_LIMIT'
    elif response['status'] not in ('NEW','PREPARED'):
        print(json.dumps(response,ensure_ascii=False));return
    payload=BASE/'payload'
    for name in ['common.py','remote_build.py','mpi_check.py']:shutil.copy2(BASE/name,payload/name)
    shutil.copy2(BASE/'build_plan.json',payload/'build_plan.json')
    argv=['rsync','-az','--no-owner','--no-group','--no-devices','--no-specials','-e',shlex.join(SSH[:-1]),str(payload)+'/', 'vast-mirheo:'+build+'/tools/']
    r=subprocess.run(argv,capture_output=True,text=True);write(BASE/'records/deploy_transfer.json',dict(argv=argv,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr));assert r.returncode==0
    response=ssh_code('''import subprocess,json,sys,shlex,os
from pathlib import Path
r=json.loads(sys.argv[1]);p=Path(r['build']);session=r['session']
check=subprocess.run(['tmux','has-session','-t',session],capture_output=True)
assert check.returncode!=0,'BUILD_SESSION_ALREADY_RUNNING'
command=shlex.join(['/usr/bin/python3','-B',str(p/'tools/remote_build.py'),str(p/'tools/build_plan.json'),str(r['attempt'])])+' > '+shlex.quote(str(p/('driver_a'+str(r['attempt'])+'.log')))+' 2>&1'
subprocess.run(['tmux','new-session','-d','-s',session,command],check=True)
print(json.dumps({'status':'BUILD_STARTED','session':session,'build_dir':str(p)}))
''',dict(build=build,session='bf-'+plan['build_id'],attempt=attempt))
    write(BASE/'records/build_launch.json',response);print(json.dumps(response))

def status(job_id=None):
    if job_id:
        remote_dir,_=result_location(job_id)
        response=ssh_code('''from pathlib import Path
import json,sys
r=json.loads(sys.argv[1]);p=Path(r['root']);s=json.loads((p/'run_status.json').read_text()) if (p/'run_status.json').is_file() else {'status':'PREPARED'}
s['results_ready']=(p/'RESULTS_READY').is_file()
if (p/'live_progress.json').exists():s['live_progress']=json.loads((p/'live_progress.json').read_text())
if s.get('solver_pid'):
 proc=Path('/proc')/str(s['solver_pid']);s['solver_pid_exists']=proc.exists()
 if proc.exists():
  try:s['solver_command_matches']=p.name in (proc/'cmdline').read_bytes().decode(errors='replace')
  except (FileNotFoundError,PermissionError):s['solver_command_matches']=False
print(json.dumps(s))
''',dict(root=remote_dir))
        write(BASE/'records'/(job_id+'-status.json'),response);print(json.dumps(response,ensure_ascii=False));return response
    plan=read(BASE/'build_plan.json');response=ssh_code('''import json,sys,os,subprocess,shutil
from pathlib import Path
r=json.loads(sys.argv[1]);p=Path(r['build']);s=json.loads((p/'build_status.json').read_text()) if (p/'build_status.json').exists() else {'status':'NOT_STARTED'}
pid=s.get('active_child_pid');s['active_child_exists']=bool(pid and Path('/proc/'+str(pid)).exists());s['workspace_free_bytes']=shutil.disk_usage('/workspace').free
print(json.dumps(s))
''',dict(build=plan['build_dir']))
    write(BASE/'records/build_status_latest.json',response);print(json.dumps(response,ensure_ascii=False))

def result_location(job_id,build=False):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{1,90}',job_id):raise ValueError('UNSAFE_JOB_ID')
    plan=read(BASE/'build_plan.json')
    if build:
        if not job_id.startswith(plan['build_id']):
            candidates=[read(p) for p in (BASE/'records').glob('build_plan_*.json') if job_id.startswith(read(p)['build_id'])]
            assert len(candidates)==1,'UNKNOWN_BUILD_ID';plan=candidates[0]
        suffix=job_id.removeprefix(plan['build_id'])
        assert job_id.startswith(plan['build_id']) and suffix in ('','-a2','-a3'),'UNKNOWN_BUILD_ID'
        remote=plan['build_dir']+'/evidence'+suffix.replace('-','_')
    else:remote='/workspace/bloodflow/cloud_runs/'+job_id
    return remote,Path('/home/lzy/projects/cloud_results')/job_id

def result_probe(remote_dir,preview=False):
    code=(BASE/'common.py').read_text()+'''
import sys
r=json.loads(sys.argv[1]);root=Path(r['root'])
allowed=root.parent==Path('/workspace/bloodflow/cloud_runs') or (root.name.startswith('evidence') and root.parent.parent==Path('/workspace/bloodflow/build'))
assert allowed and root.resolve()==root,'UNSAFE_REMOTE_RESULT_DIRECTORY'
check=validate_result(root);manifest=read(root/'results_manifest.json');ready=read(root/'RESULTS_READY')
result=dict(check=check,manifest=manifest,ready=ready,manifest_bytes=(root/'results_manifest.json').read_text(),ready_bytes=(root/'RESULTS_READY').read_text())
if r['preview']:
 users=[];inaccessible=[]
 for proc in Path('/proc').iterdir():
  if not proc.name.isdigit():continue
  try:
   for fd in (proc/'fd').iterdir():
    try:
     target=os.readlink(fd)
     if target==str(root) or target.startswith(str(root)+'/'):users.append(dict(pid=int(proc.name),path=target))
    except (FileNotFoundError,ProcessLookupError):pass
  except PermissionError:inaccessible.append(proc.name)
  except (FileNotFoundError,ProcessLookupError):pass
 result.update(open_file_users=users,inaccessible_processes=inaccessible,candidates=manifest['files'] if not users and not inaccessible else [],deleted_files=0)
print(json.dumps(result))
'''
    return ssh_code(code,dict(root=remote_dir,preview=preview))

def fetch(job_id,build=False):
    local_environment();remote_dir,dest=result_location(job_id,build)
    result=result_probe(remote_dir)
    if dest.exists():
        checked=validate_result(dest)
        assert checked['manifest_sha256']==result['check']['manifest_sha256'],'REMOTE_RESULTS_CHANGED_AFTER_ARCHIVE'
        print(json.dumps(dict(**checked,already_archived=True,local_dir=str(dest))));return checked
    parent=dest.parent;parent.mkdir(parents=True,exist_ok=True)
    partial=parent/(job_id+'.partial')
    assert parent.resolve()==parent and not partial.is_symlink(),'UNSAFE_LOCAL_RESULT_DIRECTORY'
    partial.mkdir(exist_ok=True)
    selected=[safe_relative(x['path']) for x in result['manifest']['files']]+['results_manifest.json','RESULTS_READY']
    assert len(selected)==len(set(selected))
    listing=BASE/'records'/(job_id+'-fetch-files.nul');listing.write_bytes(b''.join(x.encode()+b'\0' for x in selected))
    argv=['rsync','-rtpz','--partial','--no-links','--no-owner','--no-group','--no-devices','--no-specials','--from0','--files-from='+str(listing),
        '-e',shlex.join(SSH[:-1]),'vast-mirheo:'+remote_dir+'/',str(partial)+'/']
    r=subprocess.run(argv,capture_output=True,text=True);write(BASE/'records'/(job_id+'-fetch.json'),dict(recorded_at=now(),argv=argv,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr))
    if r.returncode:raise RuntimeError('FETCH_INCOMPLETE: same JOB_ID partial retained')
    checked=validate_result(partial)
    assert checked['manifest_sha256']==result['check']['manifest_sha256'],'RESULT_MANIFEST_CHANGED_DURING_FETCH'
    # Parse actual structured output. Failed jobs still archive their closed logs and explicit status.
    structured=[]
    for p in partial.rglob('*.json'):
        read(p);structured.append(p.relative_to(partial).as_posix())
    import csv,math
    csv_rows={}
    for p in partial.rglob('*.csv'):
        with p.open() as f:rows=list(csv.reader(f))
        csv_rows[p.relative_to(partial).as_posix()]=max(0,len(rows)-1)
    checked.update(verified_at=now(),remote_dir=remote_dir,structured_json_files=len(structured),csv_data_rows=csv_rows)
    write(partial/'LOCAL_ARCHIVE_VERIFIED.json',checked)
    partial.rename(dest)
    print(json.dumps(dict(**checked,local_dir=str(dest)),ensure_ascii=False));return checked

def verify(job_id,build=False):
    _,dest=result_location(job_id,build);result=validate_result(dest)
    print(json.dumps(result));return result

def cleanup(job_id,dry_run):
    if not dry_run:raise ValueError('DELETION_NOT_AUTHORIZED: this tool implements preview only')
    remote_dir,dest=result_location(job_id);local=validate_result(dest);remote=result_probe(remote_dir,preview=True)
    assert local['manifest_sha256']==remote['check']['manifest_sha256'],'REMOTE_CHANGED_NO_CLEANUP_CANDIDATES'
    preview=dict(job_id=job_id,dry_run=True,local_archive_verified=True,remote_hashes_unchanged=True,
        candidates=remote['candidates'],open_file_users=remote['open_file_users'],inaccessible_processes=remote['inaccessible_processes'],
        verified_files_waiting_for_process_clearance=remote['manifest']['files'] if remote['inaccessible_processes'] or remote['open_file_users'] else [],
        status='PREVIEW_BLOCKED_BY_PROCESS_VISIBILITY' if remote['inaccessible_processes'] else 'PREVIEW_BLOCKED_BY_OPEN_FILES' if remote['open_file_users'] else 'PREVIEW_READY',deleted_files=0)
    write(BASE/'records'/(job_id+'-cleanup-preview.json'),preview)
    print(json.dumps(dict(job_id=job_id,candidate_files=len(preview['candidates']),candidate_bytes=sum(x['size'] for x in preview['candidates']),deleted_files=0,record=str(BASE/'records'/(job_id+'-cleanup-preview.json')))))

def smoke(new=False,repair=True):
    if not new:raise ValueError('NEW_JOB_REQUIRED: use smoke --new; fetch/verify never rerun a solver')
    local_environment();plan=read(BASE/'build_plan.json')
    remote=ssh_code('''from pathlib import Path
import json,sys,subprocess,shutil,hashlib
r=json.loads(sys.argv[1]);b=Path(r['build']);s=json.loads((b/'build_status.json').read_text());assert s['status']=='CLOUD_BUILD_PASS','BUILD_NOT_PASSED'
e=Path(s['evidence_directory']);i=json.loads((e/'compile_identity.json').read_text());assert hashlib.sha256(Path(i['runtime_binary_path']).read_bytes()).hexdigest()==i['runtime_binary_hash']
gpu=subprocess.run(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],capture_output=True,text=True,timeout=15);assert gpu.returncode==0 and not gpu.stdout.strip(),'OTHER_GPU_PROCESS_PRESENT'
checks=[]
for f in Path('/workspace/bloodflow/build').glob('*/evidence*/build_status.json'):
 old=json.loads(f.read_text())
 for c in old['commands']:
  if c['label'].startswith('mpi_'):checks.append(dict(source=str(f),name=c['label'],elapsed_s=c['elapsed_s'],exit_code=c['exit_code']))
for f in Path('/workspace/bloodflow/build').glob('*/mpi_diagnostic*/result.json'):
 c=json.loads(f.read_text());checks.append(dict(source=str(f),name='mpi_diagnostic',elapsed_s=c['elapsed_s'],exit_code=c['exit_code']))
print(json.dumps(dict(build_status=s,identity=i,previous_checks=checks,disk_free=shutil.disk_usage('/workspace').free)))
''',dict(build=plan['build_dir']))
    assert disk_allowed(remote['disk_free'],'run')
    identity=remote['identity'];assert identity['original_source_sha256']==plan['original_source_sha256']
    job_id='cloud-smoke-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');root='/workspace/bloodflow/cloud_runs/'+job_id
    native=read(UPLOAD/'metadata/native_source_locations.json');mapping=read(UPLOAD/'metadata/path_mapping.json')
    source_mapping=[x for x in mapping['roots'] if native['latest_spec']['path'].startswith(x['local']+'/')][0]
    relative=Path(native['latest_spec']['path']).relative_to(source_mapping['local']).as_posix()
    original_spec=UPLOAD/source_mapping['staging']/relative;s=read(original_spec)
    tools_dir=BASE/'jobs'/job_id/'tools';tools_dir.mkdir(parents=True,exist_ok=False)
    for name in ['common.py','remote_smoke.py','smoke_worker.py','cloud_geometry.py','cloud_repair_worker.py','cloud_worker_adaptation.patch']:
        shutil.copy2(BASE/name,tools_dir/name)
    job_root=tools_dir.parent;tool_hashes={p.name:sha(p) for p in tools_dir.iterdir()}
    checks=remote['previous_checks'];prior=sum(x['elapsed_s'] for x in checks)
    cfg=dict(job_id=job_id,build_id=plan['build_id'],created_at=now(),original_source_path=plan['original_source_path'],resolved_cloud_path=plan['resolved_cloud_source'],
        source_hash=plan['original_source_sha256'],runtime_binary_hash=identity['runtime_binary_hash'],runtime_binary_path=identity['runtime_binary_path'],python_package=identity['python_package'],
        python=plan['python'],cloud_project_root=plan['cloud_project_root'],output_directory=root,mpi_ranks=2,compute_ranks=1,postprocess_ranks=1,gpus=1,
        child_environment=dict(HWLOC_COMPONENTS='-opencl,-gl',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1'),
        necessary_mesh=dict(original_source_path=s['mesh'],resolved_cloud_path=source_mapping['remote']+'/'+str(Path(s['mesh']).relative_to(source_mapping['local'])),source_hash=sha(UPLOAD/source_mapping['staging']/Path(s['mesh']).relative_to(source_mapping['local']))),
        parameter_source=dict(path=native['latest_spec']['path'],sha256=sha(original_spec)),liquid_parameters=s['config']['dpd'],dt=s['dt'],seed=s['seed'],liquid_domain=[4.,4.,4.],liquid_steps=100,
        run_repair_path=repair,tool_hashes=tool_hashes,limits=plan['limits'],previous_checks=checks,previous_checks_elapsed_s=prior,
        budget_ledger='/workspace/bloodflow/metadata/cloud_smoke_budget_20260910.json',authority='Current user cloud deployment and small-test request; independent of historical RBC authorization',
        physical_model_acceptance='UNCHANGED: NOT_MATCHED; no long endpoint or physical parameter search')
    cfg['cloud_plan_sha256']=object_hash(cfg)
    spec=dict(s);spec.pop('plan_sha256',None);spec['config']=json.loads(json.dumps(s['config']))
    spec.update(mesh=cfg['necessary_mesh']['resolved_cloud_path'],native_library=identity['runtime_binary_path'],native_library_sha256=identity['runtime_binary_hash'],
        cloud_plan_sha256=cfg['cloud_plan_sha256'],prep_steps=250,steps=0,sample_steps=250,continuous_observation=True,stop_after_relaxation=True,diagnostic_state=True)
    spec['config']['protocol']['relaxation_time']=250*spec['dt']
    write(job_root/'cloud_runtime.json',cfg);write(job_root/'repair_spec.json',spec);write(job_root/'compile_identity.json',identity)
    write(job_root/'worker_derivation.json',read(BASE/'records/worker_derivation.json'))
    (job_root/'cloud_build_compatibility.patch').write_text(plan['compatibility_patch'])
    write(BASE/'cloud_runtime.json',cfg)
    print(json.dumps(dict(status='CLOUD_SMOKE_PLANNED',job_id=job_id,liquid_steps=100,repair_steps=250 if repair else 0,wall_preparation_steps=4000 if repair else 0,
        cumulative_prior_checks_s=prior,limits=cfg['limits'],remote_dir=root),ensure_ascii=False),flush=True)
    ssh_code('''from pathlib import Path
import json,sys
r=json.loads(sys.argv[1]);p=Path(r['root']);assert p.parent==Path('/workspace/bloodflow/cloud_runs') and p.resolve()==p
p.mkdir(parents=True,exist_ok=False);print(json.dumps({'status':'RESERVED'}))
''',dict(root=root))
    argv=['rsync','-az','--no-owner','--no-group','--no-devices','--no-specials','-e',shlex.join(SSH[:-1]),str(job_root)+'/', 'vast-mirheo:'+root+'/']
    r=subprocess.run(argv,capture_output=True,text=True);write(BASE/'records'/(job_id+'-input_transfer.json'),dict(argv=argv,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr));assert r.returncode==0
    launch=ssh_code('''from pathlib import Path
import json,sys,subprocess,shlex
r=json.loads(sys.argv[1]);p=Path(r['root']);session='bf-'+p.name
command=shlex.join(['/usr/bin/python3','-B',str(p/'tools/remote_smoke.py'),str(p/'cloud_runtime.json')])+' > '+shlex.quote('/workspace/bloodflow/metadata/'+p.name+'-runner.log')+' 2>&1'
subprocess.run(['tmux','new-session','-d','-s',session,command],check=True);print(json.dumps({'status':'CLOUD_SMOKE_STARTED','session':session,'job_id':p.name}))
''',dict(root=root))
    write(BASE/'records'/(job_id+'-launch.json'),launch);print(json.dumps(launch),flush=True)
    # The remote tmux runner does not depend on this local polling loop staying online.
    until=time.monotonic()+210
    while time.monotonic()<until:
        result=ssh_code('''from pathlib import Path
import json,sys
p=Path(json.loads(sys.argv[1])['root']);print(json.dumps(dict(ready=(p/'RESULTS_READY').exists(),status=json.loads((p/'run_status.json').read_text()) if (p/'run_status.json').exists() else None)))
''',dict(root=root))
        if result['ready']:
            write(BASE/'records'/(job_id+'-completed.json'),result);return fetch(job_id)
        time.sleep(3)
    raise RuntimeError('REMOTE_STATUS_PENDING: use status/fetch for the same JOB_ID; do not rerun')

def rbc_full(preflight_only=False,new=False,execute=False):
    from rbc_plan import preflight,authorization
    if preflight_only:
        if new or execute:raise ValueError('PREFLIGHT_CANNOT_EXECUTE')
        return preflight()
    if not (new and execute):raise ValueError('EXPLICIT_NEW_AND_EXECUTE_REQUIRED')
    # An execution flag never invents a separate purpose/budget approval.
    frozen=read(BASE/'rbc_full/frozen_plan.json');approved=authorization(BASE,frozen)
    plan,probe=preflight()
    assert plan['plan_sha256']==frozen['plan_sha256'],'FROZEN_PLAN_CHANGED'
    assert not probe['blockers'],'PREFLIGHT_BLOCKED'
    authorization(BASE,plan)
    job='cloud-rbc-full-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root='/workspace/bloodflow/cloud_runs/'+job;local=BASE/'jobs'/job;local.mkdir(parents=True,exist_ok=False)
    (local/'tools').mkdir()
    for name,h in plan['tool_hashes'].items():
        assert sha(BASE/name)==h,'WORKER_CHANGED_AFTER_PREFLIGHT';shutil.copy2(BASE/name,local/'tools'/name)
    spec=json.loads(json.dumps(plan['spec']));spec['cloud_plan_sha256']=plan['plan_sha256']
    write(local/'full_spec.json',spec)
    mesh=[r for r in plan['python_inputs'] if r['remote']==spec['mesh']]
    assert len(mesh)==1 and sha(mesh[0]['local'])==mesh[0]['sha256'],'REFERENCE_MESH_CHANGED'
    shutil.copy2(mesh[0]['local'],local/'reference.off')
    cfg=dict(plan,job_id=job,output_directory=root,authorization=approved,spec_sha256=sha(local/'full_spec.json'),environment_before=probe['cloud']['project_venv'])
    write(local/'cloud_runtime.json',cfg);write(local/'actual_parameters.json',dict(config=spec['config'],dt=spec['dt'],counts=plan['counts'],preparation_criteria=spec['preparation_criteria']))
    write(local/'provenance.json',dict(plan_sha256=plan['plan_sha256'],identity=probe['cloud']['identity'],source_path=plan['original_source_path'],source_hash=plan['original_source_sha256'],
                                      python_inputs=plan['python_inputs'],tool_hashes=plan['tool_hashes'],worker_derivation=read(BASE/'rbc_full_worker_derivation.json'),handoff=plan['handoff'],environment=plan['environment']))
    write(local/'preflight.json',probe)
    # Reserve this authorization locally before network mutations. A lost connection cannot cause a second launch.
    receipt=BASE/'rbc_full/launch_receipt.json'
    with receipt.open('x') as f:json.dump(dict(job_id=job,plan_sha256=plan['plan_sha256'],reserved_at=now(),status='RESERVED',remote_dir=root),f,indent=2)
    ssh_code('''from pathlib import Path
import json,sys
p=Path(json.loads(sys.argv[1])['root']);assert p.parent==Path('/workspace/bloodflow/cloud_runs') and p.resolve()==p
p.mkdir(exist_ok=False);print(json.dumps({'reserved':str(p)}))
''',dict(root=root))
    argv=['rsync','-az','--no-links','--no-owner','--no-group','--no-devices','--no-specials','-e',shlex.join(SSH[:-1]),str(local)+'/', 'vast-mirheo:'+root+'/']
    r=subprocess.run(argv,capture_output=True,text=True);write(BASE/'records'/(job+'-input_transfer.json'),dict(argv=argv,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr));assert r.returncode==0
    launch=ssh_code('''from pathlib import Path
import json,sys,subprocess,shlex,hashlib
r=json.loads(sys.argv[1]);p=Path(r['root']);cfg=json.loads((p/'cloud_runtime.json').read_text())
for name,h in cfg['tool_hashes'].items():assert hashlib.sha256((p/'tools'/name).read_bytes()).hexdigest()==h
assert not (p/'STARTED').exists(),'ALREADY_STARTED'
command=shlex.join([cfg['python'],'-B',str(p/'tools/remote_rbc.py'),str(p/'cloud_runtime.json')])+' > '+shlex.quote('/workspace/bloodflow/metadata/'+p.name+'-runner.log')+' 2>&1'
subprocess.run(['tmux','new-session','-d','-s','bf-'+p.name,command],check=True)
print(json.dumps({'job_id':p.name,'status':'LAUNCH_SUBMITTED','remote_dir':str(p)}))
''',dict(root=root))
    write(BASE/'records'/(job+'-launch.json'),launch);print(json.dumps(launch),flush=True)
    until=time.monotonic()+1950
    while time.monotonic()<until:
        result=status(job)
        if result['results_ready']:
            checked=fetch(job);verify(job)
            review_rbc(job)
            return checked
        time.sleep(10)
    raise RuntimeError('REMOTE_STATUS_PENDING: status/fetch same JOB_ID; do not launch again')


def review_rbc(job):
    _,archive=result_location(job)
    destination=archive.with_name(job+'-review')
    code='import sys;sys.path.insert(0,sys.argv[1]);from rbc_report import render_review;print(render_review(sys.argv[2],sys.argv[3]))'
    r=subprocess.run(['/home/lzy/projects/mirheo_starter/.venv/bin/python','-B','-c',code,str(BASE),str(archive),str(destination)],capture_output=True,text=True)
    write(BASE/'records'/(job+'-review.json'),dict(exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr))
    if r.returncode:raise RuntimeError('REVIEW_EXPORT_FAILED: existing raw archive retained; no solver retry')
    print(r.stdout,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest='action',required=True)
    p=subs.add_parser('status');p.add_argument('job_id',nargs='?')
    p=subs.add_parser('preflight');p.add_argument('--revise-build-plan',action='store_true')
    p=subs.add_parser('deploy');p.add_argument('--resume',action='store_true',help='Explicit reviewed resume after a failed build; retains failed evidence and reuses the same work/build')
    for action in ['fetch','verify']:
        p=subs.add_parser(action);p.add_argument('job_id');p.add_argument('--build',action='store_true')
    p=subs.add_parser('cleanup');p.add_argument('job_id');p.add_argument('--dry-run',action='store_true')
    p=subs.add_parser('smoke');p.add_argument('--new',action='store_true');p.add_argument('--skip-repair',action='store_true')
    p=subs.add_parser('rbc-full');p.add_argument('--preflight-only',action='store_true');p.add_argument('--new',action='store_true');p.add_argument('--execute',action='store_true')
    p=subs.add_parser('rbc-review');p.add_argument('job_id',help='Render from an existing verified archive; never starts a solver')
    args=parser.parse_args()
    if args.action in ('fetch','verify'):{'fetch':fetch,'verify':verify}[args.action](args.job_id,args.build)
    elif args.action=='cleanup':cleanup(args.job_id,args.dry_run)
    elif args.action=='deploy':deploy(args.resume)
    elif args.action=='preflight':preflight(args.revise_build_plan)
    elif args.action=='smoke':smoke(args.new,not args.skip_repair)
    elif args.action=='rbc-full':rbc_full(args.preflight_only,args.new,args.execute)
    elif args.action=='rbc-review':review_rbc(args.job_id)
    else:status(args.job_id)
if __name__=='__main__':
    try:main()
    except Exception as e:print(json.dumps(dict(status='STOPPED',error=type(e).__name__+': '+str(e)),ensure_ascii=False));sys.exit(2)

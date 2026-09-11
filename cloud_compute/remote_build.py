"""One detached, bounded build. No source patching outside the designated work copy."""
from common import *
import subprocess,shutil,signal,sys,resource,re

def check_upload(plan):
    root=Path(plan['upload_dir']);errors=[]
    for x in plan['upload_protection']:
        p=root/safe_relative(x['path'])
        if x.get('type')=='symlink':
            if not p.is_symlink() or os.readlink(p)!=x['target']:errors.append(x['path'])
        elif p.is_symlink() or not p.is_file() or sha(p)!=x['sha256'] or p.stat().st_mode&0o111!=x['execute_bits']:errors.append(x['path'])
    actual=[]
    for parent,dirs,names in os.walk(root,followlinks=False):
        for name in dirs+names:
            p=Path(parent)/name
            if p.is_file() or p.is_symlink():actual.append(p.relative_to(root).as_posix())
    if set(actual)!={x['path'] for x in plan['upload_protection']}:errors.append('FILE_SET_CHANGED')
    if errors:raise RuntimeError('UPLOAD_SNAPSHOT_CHANGED:'+str(errors[:10]))
    return dict(status='PASS',files=len(actual),fingerprint=object_hash(plan['upload_protection']))

def readonly_environment():
    rows=[];root=Path('/venv/main')
    for base,dirs,names in os.walk(root,followlinks=False):
        for name in dirs+names:
            p=Path(base)/name;s=p.lstat()
            rows.append((str(p.relative_to(root)),s.st_size,s.st_mtime_ns,s.st_ctime_ns,stat.S_IMODE(s.st_mode),os.readlink(p) if p.is_symlink() else ''))
    return dict(path=str(root),entries=len(rows),stat_fingerprint=object_hash(sorted(rows)))

def main():
    plan=read(sys.argv[1]);b=Path(plan['build_dir']);attempt=int(sys.argv[2]) if len(sys.argv)>2 else 1
    e=b/('evidence' if attempt==1 else 'evidence_a'+str(attempt));e.mkdir(exist_ok=True)
    assert not (e/'RESULTS_READY').exists(),'COMPLETED_BUILD_EVIDENCE_IS_IMMUTABLE'
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    prior=read(b/'build_status.json') if (b/'build_status.json').exists() else {}
    charged_before=prior.get('cumulative_elapsed_s',prior.get('elapsed_s',0)) if attempt>1 else plan.get('prior_deployment_elapsed_s',0)
    t=time.monotonic();deadline=t+plan['limits']['build_seconds']-charged_before;records=[]
    env=dict(os.environ,HWLOC_COMPONENTS='-opencl,-gl',PYTHONDONTWRITEBYTECODE='1',PIP_DISABLE_PIP_VERSION_CHECK='1',
        CUDA_VISIBLE_DEVICES='0',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',TMPDIR=str(b/'tmp'),
        PATH='/workspace/bloodflow/.venv/bin:/usr/local/cuda/bin:/usr/bin:/bin',LD_LIBRARY_PATH='/usr/local/cuda/lib64')
    env.pop('PYTHONPATH',None);(b/'tmp').mkdir(exist_ok=True)
    state=dict(build_id=plan['build_id'],attempt=attempt,status='RUNNING',pid=os.getpid(),started_at=now(),commands=records,build_plan_sha256=sha(sys.argv[1]),
        evidence_directory=str(e),environment_overrides={'HWLOC_COMPONENTS':'-opencl,-gl'},startup_fix='Disable GL topology plugin as well as OpenCL; verbose first failure showed GL enabled before hang; no solver changes')
    write(b/'build_status.json',state)
    def run(label,argv,extra_env=None,limit=None):
        if not disk_allowed(shutil.disk_usage('/workspace').free,'run'):raise RuntimeError('DISK_LIMIT')
        start=time.monotonic();log=e/(label+'.log');timeout=min(deadline-start,limit or 1800)
        if timeout<=0:raise TimeoutError('BUILD_TIME_LIMIT')
        with log.open('wb') as f:
            child=subprocess.Popen(argv,stdout=f,stderr=subprocess.STDOUT,env={**env,**(extra_env or {})},start_new_session=True,cwd=b)
            state['active_child_pid']=child.pid;state['active_command']=label;write(b/'build_status.json',state)
            stop=None
            while child.poll() is None:
                if time.monotonic()-start>=timeout:stop='TIMEOUT'
                elif not disk_allowed(shutil.disk_usage('/workspace').free,'run'):stop='DISK_LIMIT'
                if stop:
                    os.killpg(child.pid,signal.SIGTERM)
                    try:child.wait(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                    break
                time.sleep(1)
        record=dict(label=label,argv=argv,exit_code=child.returncode,elapsed_s=time.monotonic()-start,stop_reason=stop,log=str(log))
        records.append(record);state['active_child_pid']=None;write(b/'build_status.json',state)
        if child.returncode or stop:raise RuntimeError(label+':'+str(child.returncode)+':'+str(stop))
        return log
    final='FAILED'
    try:
        assert disk_allowed(shutil.disk_usage('/workspace').free,'build',plan['estimated_build_growth_bytes']),'INSUFFICIENT_BUILD_SPACE'
        write(e/'upload_before.json',check_upload(plan));write(e/'image_environment_before.json',readonly_environment())
        tool_before={p:sha(p) for p in plan['protected_tool_binaries']};write(e/'tool_binaries_before.json',tool_before)
        python=plan['python'];run('packages_before',[python,'-B','-m','pip','list','--format=json'],limit=20)
        source=Path(plan['resolved_cloud_source']);work=Path(plan['work_source'])
        if not work.exists() and plan.get('work_source_reuse_from'):
            old=Path(plan['work_source_reuse_from']);assert old.parent.parent==Path('/workspace/bloodflow/work') and old.name=='source'
            assert old.is_dir() and not old.is_symlink();work.parent.mkdir(parents=True,exist_ok=True)
            old.rename(work)
            write(e/'work_source_relocation.json',dict(from_path=str(old),to_path=str(work),reason='Reuse the single writable source copy after the archived CMake architecture-parser failure; no input/result deletion'))
            for item in plan['compatibility_edits']:
                p=work/item['path'];h=sha(p)
                if h==item['original_sha256']:p.write_text(item['replacement'])
                else:assert h==item['replacement_sha256'],'UNRECOGNIZED_REUSED_SOURCE_EDIT'
        if not work.exists():
            shutil.copytree(source,work,symlinks=True)
            for item in plan['compatibility_edits']:
                p=work/safe_relative(item['path']);assert sha(p)==item['original_sha256'],'COMPATIBILITY_PATCH_SOURCE_MISMATCH'
                p.write_text(item['replacement']);assert sha(p)==item['replacement_sha256']
        else:
            assert attempt>1 or plan.get('work_source_reuse_from'),'EXISTING_WORK_SOURCE_REQUIRES_EXPLICIT_REVIEWED_RESUME'
            for item in plan['compatibility_edits']:assert sha(work/item['path'])==item['replacement_sha256']
        write(e/'cloud_compatibility_edits.json',plan['compatibility_edits'])
        (e/'cloud_build_compatibility.patch').write_text(plan['compatibility_patch'])
        for relative,original_sha in plan['protected_physics_files'].items():assert sha(work/relative)==original_sha,'PHYSICS_SOURCE_CHANGED'
        run('mpi_hostname',['mpirun','--allow-run-as-root','--bind-to','none','-np','2','hostname'],limit=10)
        run('mpi_communication',['mpirun','--allow-run-as-root','--bind-to','none','-np','2',python,'-B',str(b/'tools/mpi_check.py'),plan['cloud_project_root']],limit=20)
        # Only missing packages are installed; already satisfied versions remain untouched.
        installed={x['name'].lower():x['version'] for x in json.loads((e/'packages_before.log').read_text())}
        missing=[s for s in plan['python_dependencies'] if installed.get(s.split('==')[0].lower())!=s.split('==')[1]]
        write(e/'dependency_changes_requested.json',dict(requested=missing,existing=installed))
        if missing:run('install_dependencies',[python,'-B','-m','pip','--isolated','install','--no-cache-dir','--index-url','https://pypi.org/simple','--report',str(e/'dependencies_install_report.json'),*missing],limit=180)
        run('configure',plan['configure_argv'])
        cache=(b/'CMakeCache.txt').read_text();flags=re.search(r'^MIR_CUDA_ARCH_FLAGS:INTERNAL=(.*)$',cache,re.M)
        assert flags and 'compute_120' in flags[1] and 'code=sm_120' in flags[1] and 'sm_89' not in flags[1],'WRONG_MIRHEO_CUDA_ARCHITECTURE'
        assert 'HDF5_C_LIBRARY_hdf5:FILEPATH=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so' in cache,'PARALLEL_HDF5_NOT_SELECTED'
        run('compile',['cmake','--build',str(b),'--parallel','2','--verbose'])
        libraries=list((b/'src/mirheo/bindings').glob('libmirheo.cpython*.so'));assert len(libraries)==1
        library=libraries[0];commands=read(b/'compile_commands.json');cuda=[x for x in commands if x['file'].endswith('.cu')]
        assert cuda and all('compute_120' in x['command'] and 'sm_120' in x['command'] for x in cuda),'COMPILE_COMMAND_TARGET_NOT_120'
        targets=set(re.findall(r'(?:sm|compute)_(\d+)', '\n'.join(x['command'] for x in cuda)));assert targets=={'120'},targets
        cuobjdump=run('cuobjdump',['/usr/local/cuda/bin/cuobjdump','--list-elf',str(library)],limit=45)
        binary_targets=set(re.findall(r'sm_(\d+)',cuobjdump.read_text()));assert binary_targets=={'120'},binary_targets
        linklog=run('ldd',['ldd',str(library)],limit=20);assert 'not found' not in linklog.read_text(),'DYNAMIC_DEPENDENCY_MISSING'
        run('install_mirheo',[python,'-B','-m','pip','--isolated','install','--no-deps','--no-build-isolation','--no-cache-dir','--report',str(e/'mirheo_install_report.json'),str(work)],extra_env={'MIRHEO_PREBUILT_LIBRARY':str(library)},limit=120)
        site=run('python_site',[python,'-B','-c','import sysconfig;print(sysconfig.get_path("platlib"))'],limit=20).read_text().strip()
        site=Path(site);assert site.is_relative_to(Path('/workspace/bloodflow/.venv'))
        installed_library=site/library.name;assert sha(installed_library)==sha(library),'INSTALLED_BINARY_HASH_MISMATCH'
        run('packages_after',[python,'-B','-m','pip','list','--format=json'],limit=20)
        write(e/'compile_identity.json',dict(build_id=plan['build_id'],original_source_sha256=plan['original_source_sha256'],
            runtime_binary_path=str(installed_library),runtime_binary_hash=sha(installed_library),build_binary_path=str(library),
            python_package=str(site/'mirheo'),target_architectures=sorted(binary_targets),cuda_compile_units=len(cuda),
            precision_bits=32,membrane_double=False,rod_double=False,build_type='Release',fast_math=True,
            source_revision=plan['source_revision'],compatibility_patch_sha256=hashlib.sha256(plan['compatibility_patch'].encode()).hexdigest(),
            archive_patch_reapplied=False,physics_sources_unchanged=True))
        shutil.copy2(b/'CMakeCache.txt',e/'CMakeCache.txt');shutil.copy2(b/'compile_commands.json',e/'compile_commands.json')
        for relative,original_sha in plan['protected_physics_files'].items():assert sha(work/relative)==original_sha,'PHYSICS_SOURCE_CHANGED'
        final='CLOUD_BUILD_PASS'
    except Exception as ex:
        state['error']=type(ex).__name__+': '+str(ex)
    finally:
        try:
            write(e/'upload_after.json',check_upload(plan))
            image_after=readonly_environment();write(e/'image_environment_after.json',image_after)
            if (e/'image_environment_before.json').exists():assert image_after==read(e/'image_environment_before.json'),'IMAGE_ENVIRONMENT_CHANGED'
            if (e/'tool_binaries_before.json').exists():assert {p:sha(p) for p in plan['protected_tool_binaries']}==read(e/'tool_binaries_before.json'),'TOOLCHAIN_CHANGED'
        except Exception as ex:final='FAILED';state['protection_error']=str(ex)
        disk=shutil.disk_usage('/workspace');write(e/'disk_after.json',dict(free=disk.free,used=disk.used,dirs={p:subprocess.run(['du','-sb',p],capture_output=True,text=True).stdout for p in [str(b),plan['work_source'],'/workspace/bloodflow/.venv']}))
        state.update(status=final,finished_at=now(),elapsed_s=time.monotonic()-t,cumulative_elapsed_s=charged_before+time.monotonic()-t,active_child_pid=None,commands=records)
        write(b/'build_status.json',state);write(e/'build_status.json',state)
        seal(e,plan['build_id']+('' if attempt==1 else '-a'+str(attempt)),'COMPLETED' if final=='CLOUD_BUILD_PASS' else 'FAILED')
    return 0 if final=='CLOUD_BUILD_PASS' else 2

if __name__=='__main__':sys.exit(main())

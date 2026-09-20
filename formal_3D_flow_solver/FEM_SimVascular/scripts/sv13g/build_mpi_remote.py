"""The single authorized fallback: isolated, same-version Open MPI."""
import hashlib,json,os,subprocess,time
from pathlib import Path
from mpi_remote import BASE,R,L,B,ENV,run,okay,write
assert json.loads((R/'system_mpi_tests.json').read_text())['fallback_needed']
assert not (R/'mpi_fallback_build.json').exists(), 'Only one fallback build is authorized'
src=BASE/'external/openmpi-4.1.6';prefix=BASE/'external/gpu_mpi'
manifest={'status':'BUILDING','version':'4.1.6','family':'Open MPI','prefix':str(prefix),'source':json.loads((R/'mpi_fallback_source.json').read_text()),'environment':ENV,'steps':[],'fallback_count':1}
write('mpi_fallback_build',manifest)
def build(args,name):
    start=time.monotonic()
    with (L/(name+'.log')).open('w') as f:
        p=subprocess.run(args,cwd=src,env=ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
    record={'command':args,'cwd':str(src),'exit_code':p.returncode,'wall_time_s':time.monotonic()-start,'log':'logs/'+name+'.log'}
    manifest['steps'].append(record);write('mpi_fallback_build',manifest);print(json.dumps(record),flush=True)
    if p.returncode:
        manifest.update(status='FAIL',reason='MPI_RUNTIME_UNUSABLE');write('mpi_fallback_build',manifest);raise SystemExit(1)
build(['./configure','--prefix='+str(prefix),'--disable-mpi-fortran','--disable-oshmem','--with-hwloc=internal','--with-libevent=internal','--with-pmix=internal','--without-cuda','--without-ucx','--without-ofi','--without-verbs'],'mpi_fallback_configure')
build(['make','-j8'],'mpi_fallback_make')
build(['make','install'],'mpi_fallback_install')
ENV.update(PATH=str(prefix/'bin')+':/usr/bin:/bin:/usr/local/cuda/bin',LD_LIBRARY_PATH=str(prefix/'lib'))
cc=run('local_mpi_compile',[str(prefix/'bin/mpicc'),str(B/'mpi_hello.c'),'-o',str(B/'mpi_hello_local')],limit=30)
assert okay(cc)
singleton=run('local_mpi_singleton',[str(B/'mpi_hello_local')])
recipe=[str(prefix/'bin/mpiexec')]+(['--allow-run-as-root'] if os.getuid()==0 else [])
repeats=[run('local_rank1_'+str(i),recipe+['-n','1',str(B/'mpi_hello_local')]) for i in range(1,6)]
rank2=run('local_rank2',recipe+['-n','2',str(B/'mpi_hello_local')])
passed=okay(singleton) and singleton['stdout'].strip()=='rank=0 size=1' and all(okay(t) and t['stdout'].strip()=='rank=0 size=1' for t in repeats) and okay(rank2) and set(rank2['stdout'].splitlines())=={'rank=0 size=2','rank=1 size=2'}
write('mpi_hard_gate',{'status':'PASS' if passed else 'FAIL','singleton':singleton,'rank1':repeats,'rank2':rank2,'command':recipe,'environment':ENV})
libs={p:run('local_'+p+'_ldd',['ldd',str(prefix/'bin'/p)]) for p in ('mpiexec','mpicc','mpicxx')}
binary=run('local_hello_ldd',['ldd',str(B/'mpi_hello_local')])
info=run('local_ompi_info',[str(prefix/'bin/ompi_info'),'--all'])
manifest.update(status='PASS' if passed else 'FAIL',reason=None if passed else 'MPI_RUNTIME_UNUSABLE',linked_libraries=libs,hello_ldd=binary,binary_sha256=hashlib.sha256((B/'mpi_hello_local').read_bytes()).hexdigest(),compiler=run('local_cc_version',['gcc','--version']),runtime=run('local_runtime_version',[str(prefix/'bin/mpiexec'),'--version']))
write('mpi_fallback_build',manifest)
if not passed:raise SystemExit(1)
script=BASE/'scripts/run_gpu_mpi.sh';script.parent.mkdir(exist_ok=True)
script.write_text('#!/usr/bin/env bash\nset -euo pipefail\n# Validated Open MPI 4.1.6; no scientific parameters.\nexport PATH="'+str(prefix/'bin')+':/usr/bin:/bin:/usr/local/cuda/bin"\nexport LD_LIBRARY_PATH="'+str(prefix/'lib')+'${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\nroot_args=()\nif [[ "$(id -u)" == 0 ]]; then root_args=(--allow-run-as-root); fi\nexec "'+str(prefix/'bin/mpiexec')+'" "${root_args[@]}" "$@"\n')
script.chmod(0o755)
wrapper=run('verified_wrapper',[str(script),'-n','1',str(B/'mpi_hello_local')]);assert okay(wrapper)
write('mpi_resolution',{'status':'PASS','classification':'PROJECT_LOCAL_OPENMPI_RESOLVED','implementation':'Open MPI','version':'4.1.6','prefix':str(prefix),'root':os.getuid()==0,'root_override_explicit':True,'system_stack_mismatch':False,'fallback_count':1,'working_launcher':str(script),'wrapper_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'rank1_passed':5,'rank2':'PASS','singleton':'PASS','system_issue':'strace observes hwloc GL/OpenCL plugin loading and blocking X11 socket read during topology discovery; exact plugin attribution is inferred, not a captured call stack','system_MPI_modified':False,'source_family_unchanged':True})

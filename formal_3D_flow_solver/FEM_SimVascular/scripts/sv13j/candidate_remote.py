"""One authorized CUDA candidate, through original PETSc make. Stop after winner."""
import hashlib, json, os, pwd, re, shutil, sys, tarfile, time, urllib.request
from pathlib import Path
from runner_remote import *
assert load('gcc12_environment')['status']=='PASS'
key=sys.argv[1]; order=['cuda123','cuda122','cuda121']; assert key in order
version={'cuda123':'12.3.2','cuda122':'12.2.2','cuda121':'12.1.1'}[key]
matrix=load('compatibility_matrix') if (R/'compatibility_matrix.json').exists() else {'winner':None,'candidates':[]}
assert matrix['winner'] is None,'STOP_ON_FIRST_MAKE_SUCCESS'
assert [d['key'] for d in matrix['candidates']]==order[:order.index(key)],'CANDIDATE_ORDER_VIOLATION'
assert all(d['status']=='FAIL' for d in matrix['candidates'])
d={'key':key,'CUDA_version':version,'status':'RUNNING','nvcc_version':None,'Thrust_version':None,
   'host_compiler':load('gcc12_environment')['compilers']['g++-12']['number'],
   'kernel':'NOT_RUN','MPI':'NOT_RUN','configure':'NOT_RUN','make':'NOT_RUN','first_error':None,'classification':None}
matrix['candidates'].append(d)
def save():write('compatibility_matrix',matrix)
def fail(reason,classification):
    d.update(status='FAIL',first_error=reason,classification=classification);save();print(json.dumps(d),flush=True);raise SystemExit(1)
save()
selection=json.loads((BASE/'configs'/f'{key}_selection.json').read_text())
assert selection['version']==version and selection['url'].startswith(f'https://developer.download.nvidia.com/compute/cuda/{version}/')
E=BASE/'external/compat_cuda';E.mkdir(exist_ok=True)
archive=E/selection['url'].rsplit('/',1)[1];assert not archive.exists()
start=time.monotonic();md5=hashlib.md5();sha=hashlib.sha256();size=0;notice=512*1024*1024
with urllib.request.urlopen(selection['url'],timeout=60) as src,archive.with_suffix('.part').open('xb') as out:
    length=int(src.headers['Content-Length'])
    for block in iter(lambda:src.read(8*1024*1024),b''):
        out.write(block);sha.update(block);md5.update(block);size+=len(block)
        if size>=notice:print(f'{key}: downloaded {size/1024**3:.2f}/{length/1024**3:.2f} GiB',flush=True);notice+=512*1024*1024
assert size==length and md5.hexdigest()==selection['official_md5'],'OFFICIAL_CHECKSUM_MISMATCH'
archive.with_suffix('.part').rename(archive)
write(key+'_source_manifest',dict(selection,archive=str(archive),bytes=size,actual_md5=md5.hexdigest(),sha256=sha.hexdigest(),
    frozen_before_installer_execution=True,status='PASS',download_s=time.monotonic()-start,
    checksum_provenance='NVIDIA-published MD5 verified; locally computed SHA256 frozen, not an NVIDIA-published SHA256'))
prefix=E/('cuda-'+version);tmp=BASE/('tmp_'+key);assert not prefix.exists() and not tmp.exists()
uid=pwd.getpwnam('nobody');parentstat=E.stat()
for p in (prefix,tmp):p.mkdir();os.chown(p,uid.pw_uid,uid.pw_gid)
hostbin=E/'host-bin';hostbin.mkdir(exist_ok=True)
for n in ('gcc','g++'):
    if not (hostbin/n).exists():(hostbin/n).symlink_to('/usr/bin/'+n+'-12')
cmd=['/usr/bin/setpriv','--reuid='+str(uid.pw_uid),'--regid='+str(uid.pw_gid),'--clear-groups','--no-new-privs',
     'sh',str(archive),'--silent','--toolkit','--toolkitpath='+str(prefix),'--defaultroot='+str(prefix),
     '--no-man-page','--tmpdir='+str(tmp)]
try:
    os.chown(E,uid.pw_uid,uid.pw_gid)
    install=run(cmd,key+'_install',cwd=tmp,timeout=900,extra_env={'PATH':str(hostbin)+':/usr/bin:/bin','CC':'/usr/bin/gcc-12','CXX':'/usr/bin/g++-12'})
finally:
    os.chown(E,parentstat.st_uid,parentstat.st_gid)
    for base,dirs,files in os.walk(prefix,followlinks=False):
        os.chown(base,0,0)
        for n in dirs+files:os.chown(Path(base)/n,0,0,follow_symlinks=False)
log=Path('/tmp/cuda-installer.log')
if log.exists() and str(prefix) in log.read_text(errors='replace'):shutil.copyfile(log,L/(key+'_installer_detail.log'))
write(key+'_installation',dict(install,prefix=str(prefix),toolkit_only=True,driver_install=False,source_sha256=sha.hexdigest(),status='PASS' if okay(install) else 'FAIL'))
if not okay(install):fail('CUDA_'+version+'_INSTALL_FAIL','OTHER')
mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text());mpiprefix=Path(mpi['prefix'])
wrapper=BASE/'scripts'/('use_'+key+'_env.sh')
wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\n# Scoped candidate environment; no system defaults edited.\nexport CUDA_HOME="'+str(prefix)+'"\nexport PATH="'+str(prefix/'bin')+':'+str(hostbin)+':'+str(mpiprefix/'bin')+':/usr/bin:/bin${PATH:+:$PATH}"\nexport LD_LIBRARY_PATH="'+str(prefix/'lib64')+':'+str(mpiprefix/'lib')+'${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\nexport CC=/usr/bin/gcc-12\nexport CXX=/usr/bin/g++-12\nif [[ "$#" -gt 0 ]]; then exec "$@"; fi\n')
wrapper.chmod(0o755)
env12={'OMPI_CC':'/usr/bin/gcc-12','OMPI_CXX':'/usr/bin/g++-12','PETSC_OPTIONS':'-use_gpu_aware_mpi 0'}
def execute(args,name,cwd=BASE,timeout=60):return run(args,key+'_'+name,cwd=cwd,timeout=timeout,cuda=wrapper,extra_env=env12)
nvcc=execute([prefix/'bin/nvcc','--version'],'nvcc_version');assert okay(nvcc)
assert 'release '+'.'.join(version.split('.')[:2]) in text(nvcc)
headers={}
for label,relative in [('thrust','include/thrust/version.h'),('cub','include/cub/version.cuh'),('libcudacxx','include/cuda/std/detail/libcxx/include/__config'),('cccl','include/cuda/version')]:
    p=prefix/relative
    if p.exists():
        dest=R/(key+'_headers');dest.mkdir(exist_ok=True)
        shutil.copyfile(p,dest/(label+'.txt'))
        headers[label]={'path':str(p),'sha256':digest(p),'version_macros':[line for line in p.read_text(errors='replace').splitlines() if re.match(r'^\s*#\s*define\s+\w*VERSION\w*\s+',line)]}
    else:headers[label]={'path':str(p),'exists':False,'version_macros':None}
thrust=next(x for x in headers['thrust']['version_macros'] if re.match(r'^\s*#\s*define THRUST_VERSION\s+\d+',x))
d['Thrust_version']=int(thrust.split()[-1]);d['nvcc_version']=text(nvcc)
compiler_show=[]
for n in ('mpicc','mpicxx'):
    x=execute([mpiprefix/'bin'/n,'--showme:command'],n+'_host_compiler');x['stdout']=text(x);compiler_show.append(x)
assert 'gcc-12' in compiler_show[0]['stdout'] and 'g++-12' in compiler_show[1]['stdout']
write(key+'_versions',{'status':'PASS','CUDA_version':version,'prefix':str(prefix),'nvcc_path':str(prefix/'bin/nvcc'),
    'nvcc_version':text(nvcc),'headers':headers,'THRUST_VERSION':d['Thrust_version'],'host_compiler':load('gcc12_environment'),
    'MPI_host_compiler_override':env12,'MPI_compiler_show':compiler_show,'wrapper':str(wrapper),'wrapper_sha256':digest(wrapper),
    'version_manifest':json.loads((prefix/'version.json').read_text())})
binary=B/(key+'_kernel_smoke')
compile=execute([prefix/'bin/nvcc','-ccbin','/usr/bin/g++-12','-std=c++17','-arch=sm_89','--cudart','shared',B/'cuda_kernel_smoke.cu','-o',binary],'kernel_compile',timeout=120)
runs=[]
if okay(compile):
    for i in range(1,4):
        x=execute([binary],'kernel_run_'+str(i));x['stdout']=text(x);runs.append(x)
success=len(runs)==3 and all(okay(x) and 'correct=1 last_error=0' in x['stdout'] for x in runs)
d['kernel']='PASS' if success else 'FAIL'
write(key+'_runtime',{'status':d['kernel'],'compile':compile,'runs':runs,'binary':str(binary) if success else None,'binary_sha256':digest(binary) if success else None,'source_sha256':digest(B/'cuda_kernel_smoke.cu')})
save()
if not success:fail('CUDA_'+version+'_RUNTIME_FAIL','HOST_COMPILER' if '__GNUC__' in text(compile) else 'OTHER')
assert digest(mpi['working_launcher'])==mpi['wrapper_sha256']
artifacts=json.loads((BASE/'configs/mpi_native_artifacts.json').read_text())
assert all(digest(a['realpath'])==a['sha256'] for a in artifacts)
hello=next(a['realpath'] for a in artifacts if a['path'].endswith('mpi_hello_local'))
rank1=[]
for i in range(1,6):
    x=execute([mpi['working_launcher'],'-n','1',wrapper,hello],'mpi_rank1_'+str(i),timeout=10);x['stdout']=text(x);rank1.append(x)
rank2=execute([mpi['working_launcher'],'-n','2',wrapper,hello],'mpi_rank2',timeout=10);rank2['stdout']=text(rank2)
mpi_pass=all(okay(x) and x['stdout'].strip()=='rank=0 size=1' for x in rank1) and okay(rank2) and set(rank2['stdout'].splitlines())=={'rank=0 size=2','rank=1 size=2'}
d['MPI']='PASS' if mpi_pass else 'FAIL';write(key+'_MPI',{'status':d['MPI'],'rank1':rank1,'rank2':rank2,'MPI_wrapper_unchanged':True,'MPI_rebuilt':False});save()
if not mpi_pass:fail('MPI_RUNTIME_UNUSABLE','OTHER')
source=E/('petsc-3.19.6-'+key);petscprefix=E/('petsc-'+key);arch='arch-sv13j-'+key
archive=BASE.parent/'sv1_3/petsc-3.19.6.tar.gz'
reference=json.loads((BASE/'configs/reference_manifest.json').read_text())
assert digest(archive)==reference['PETSc_archive_sha256'] and not source.exists()
original={}
with tarfile.open(archive) as tar:
    for member in tar:
        if member.isfile():original[str(Path(member.name).relative_to('petsc-3.19.6'))]=hashlib.sha256(tar.extractfile(member).read()).hexdigest()
    # Candidate-specific extraction container prevents collisions and object reuse.
    dest=E/(key+'_source_extract');dest.mkdir();tar.extractall(dest,filter='data')
(dest/'petsc-3.19.6').rename(source)
assert all(digest(source/n)==h for n,h in original.items())
write(key+'_source_integrity_before',{'status':'PASS','source_files':len(original),'archive_sha256':digest(archive),'source':str(source),'original_files':original,'modifications':[],'fresh_source':True,'PETSC_ARCH':arch,'previous_objects_reused':False})
libs=E/'system-lib';libs.mkdir(exist_ok=True)
for name in ('blas','lapack'):
    p=libs/('lib'+name+'.so')
    if not p.exists():p.symlink_to('/usr/lib/x86_64-linux-gnu/lib'+name+'.so.3')
cuda_libs=[prefix/'lib64'/('lib'+n+'.so') for n in ('cudart','cufft','cublas','cusparse','cusolver','curand','nvToolsExt')]
assert all(p.exists() for p in cuda_libs)
configure=['/usr/bin/python3','-B','configure','--prefix='+str(petscprefix),'PETSC_ARCH='+arch,
    '--with-cc='+str(mpiprefix/'bin/mpicc'),'--with-cxx='+str(mpiprefix/'bin/mpicxx'),'--with-mpiexec='+mpi['working_launcher'],
    '--with-fc=0','--with-fortran-bindings=0','--with-debugging=0','--with-shared-libraries=1',
    '--with-scalar-type=real','--with-precision=double','--with-64-bit-indices=0',
    '--with-blaslapack-lib=['+str(libs/'liblapack.so')+','+str(libs/'libblas.so')+']','--with-x=0',
    '--with-cuda=1','--with-cudac='+str(prefix/'bin/nvcc'),'--with-cuda-include='+str(prefix/'include'),
    '--with-cuda-lib=['+','.join(map(str,cuda_libs))+']','--with-cuda-arch=89',
    '--with-cxx-dialect=C++17','--with-cuda-dialect=C++17','--with-hip=0','--with-opencl=0']
cfg=execute(configure,'petsc_configure',cwd=source,timeout=3600);d['configure']='PASS' if okay(cfg) else 'FAIL'
write('petsc_'+key+'_configure',dict(cfg,status=d['configure'],source=str(source),PETSC_ARCH=arch,compiler_overrides=env12));save()
evidence=R/(key+'_build_evidence');evidence.mkdir()
if not okay(cfg):fail('PETSC_CONFIGURE_FAIL','OTHER')
variables=(source/arch/'lib/petsc/conf/petscvariables').read_text()
for name in ('CXX_FLAGS','CUDAC_FLAGS'):
    line=next(x for x in variables.splitlines() if x.startswith(name+' ='))
    assert '-std=c++17' in line and '-std=c++20' not in line
assert '#define PETSC_HAVE_CUDA 1' in (source/arch/'include/petscconf.h').read_text()
make=execute(['make','-j8','V=1','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'petsc_make',cwd=source,timeout=3600)
d['make']='PASS' if okay(make) else 'FAIL'
changed=[n for n,h in original.items() if not (source/n).is_file() or digest(source/n)!=h]
write(key+'_source_integrity_after',{'status':'FAIL' if changed else 'PASS','source_files':len(original),'modifications':changed,'archive_sha256':digest(archive)})
assert not changed,'PETSC_SOURCE_MODIFIED'
for p in (source/arch/'lib/petsc/conf/configure.log',source/arch/'lib/petsc/conf/make.log',source/arch/'lib/petsc/conf/petscvariables',source/arch/'include/petscconf.h'):
    if p.exists():shutil.copyfile(p,evidence/p.name)
errors=[{'line':i,'text':line} for i,line in enumerate(text(make).splitlines(),1) if 'error:' in line or 'undefined reference' in line or 'fatal error' in line]
unique=list(dict.fromkeys(re.sub(r'^.*?:\d+:\d+:','',x['text']) for x in errors))
full_commands=[line for line in text(make).splitlines() if str(prefix/'bin/nvcc') in line and ' -c ' in line]
manifest={'status':d['make'],'version':'3.19.6','source':str(source),'prefix':str(petscprefix),'PETSC_ARCH':arch,
    'source_archive_sha256':digest(archive),'source_unmodified':not changed,'source_files_verified':len(original),
    'CUDA_version':version,'CUDA_prefix':str(prefix),'Thrust_version':d['Thrust_version'],'host_compiler':d['host_compiler'],
    'MPI_prefix':str(mpiprefix),'configure':cfg,'make':make,'configure_command':configure,'errors':errors,
    'unique_errors':unique,'compiler_commands':full_commands,'library_sha256':digest(source/arch/'lib/libpetsc.so') if okay(make) else None}
write('petsc_'+key+'_build',manifest)
if okay(make):
    d.update(status='PASS',classification='COMPATIBLE_BUILD');matrix['winner']=key;save()
    write('compatibility_winner',dict(manifest,key=key,candidate_wrapper=str(wrapper),frozen_on_first_make_pass=True))
    print('COMPATIBILITY WINNER '+version+'; matrix search stopped',flush=True)
else:
    first=errors[0]['text'] if errors else text(make)[-2000:]
    category='THRUST_API' if any('thrust' in x.lower() for x in unique) else 'CUDA_API' if 'cuda' in first.lower() else 'LINKER' if 'undefined reference' in first else 'OTHER'
    if 'has no member "get"' in first and key=='cuda123':d['signature']='CUDA12.3_TUPLE_API_STILL_INCOMPATIBLE'
    fail(first,category)

"""Clean PETSc CUDA build, strictly gated on verified native MPI."""
import hashlib,json,os,shutil,subprocess,tarfile,time
from pathlib import Path
from mpi_remote import BASE,R,L,B,ENV,run,okay,write
resolution=json.loads((R/'mpi_resolution.json').read_text())
gate=json.loads((R/'mpi_hard_gate.json').read_text())
assert resolution['status']==gate['status']=='PASS'
assert len(gate['rank1'])==5 and all(okay(x) for x in gate['rank1']) and okay(gate['rank2'])
assert not (R/'cuda_build.json').exists(),'Do not reuse a failed build directory'
prior=json.loads((BASE/'configs/prior_cuda_build.json').read_text())
mpi=Path(resolution['prefix']);wrapper=resolution['working_launcher']
ENV.update(PATH=str(mpi/'bin')+':/usr/bin:/bin:/usr/local/cuda/bin',LD_LIBRARY_PATH=str(mpi/'lib'))
source=BASE/'external/petsc-3.19.6-cuda-clean';prefix=BASE/'external/petsc-cuda';arch='arch-sv13g-cuda'
archive=BASE.parent/'sv1_3/petsc-3.19.6.tar.gz'
sha=hashlib.sha256(archive.read_bytes()).hexdigest();assert sha==prior['source_archive_sha256']
assert not source.exists()
with tarfile.open(archive) as f:f.extractall(BASE/'external',filter='data')
(BASE/'external/petsc-3.19.6').rename(source)
lib=BASE/'external/system-lib';lib.mkdir(exist_ok=True)
for name in ('blas','lapack'):(lib/('lib'+name+'.so')).symlink_to('/usr/lib/x86_64-linux-gnu/lib'+name+'.so.3')
manifest={'status':'BUILDING','version':prior['version'],'source_archive_sha256':sha,'source':str(source),'prefix':str(prefix),'PETSC_DIR':str(source),'PETSC_ARCH':arch,'MPI_prefix':str(mpi),'MPI_verified_before_configure':True,'environment':ENV,'steps':[],'scientific_source_patched':False,'toolkit_version':json.loads((R/'cuda_version.json').read_text())['stdout']}
write('cuda_build',manifest)
# Ordinary C and MPI programs must run both directly and via verified launcher.
(B/'plain.c').write_text('int main(void) { return 0; }\n')
pre=[]
for kind in ('plain','mpi_hello'):
    exe=B/(kind+'_preconfigure')
    pre.append(run(kind+'_precompile',[str(mpi/'bin/mpicc'),str(B/(kind+'.c')),'-o',str(exe)],limit=30))
    assert okay(pre[-1])
    pre.append(run(kind+'_predirect',[str(exe)]))
    pre.append(run(kind+'_prempiexec',[wrapper,'-n','1',str(exe)]))
assert all(okay(x) for x in pre)
write('preconfigure_smoke',{'status':'PASS','tests':pre})
def step(args,name):
    start=time.monotonic();log=L/(name+'.log')
    print('START '+name,flush=True)
    with log.open('x') as f:p=subprocess.run(args,cwd=source,env=ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=3600)
    record={'name':name,'command':args,'cwd':str(source),'exit_code':p.returncode,'wall_time_s':time.monotonic()-start,'log':'logs/'+log.name,'sha256':hashlib.sha256(log.read_bytes()).hexdigest()}
    manifest['steps'].append(record);write('cuda_build',manifest)
    if name=='petsc_configure':
        if (source/'configure.log').exists():shutil.copyfile(source/'configure.log',L/'petsc_configure_detail.log')
        write('cuda_configure',dict(record,status='PASS' if p.returncode==0 else 'FAIL',version=prior['version'],MPI_verified_before_configure=True))
    print('END '+name+' exit='+str(p.returncode),flush=True)
    if p.returncode:
        manifest.update(status='FAIL',failed_step=name,reason='PENDING_ERROR_CLASSIFICATION');write('cuda_build',manifest)
        print(log.read_text(errors='replace')[-6000:],flush=True);raise SystemExit(p.returncode)
    return record
step(['/usr/bin/python3','-B','configure','--help'],'petsc_configure_help')
help=(L/'petsc_configure_help.log').read_text();assert '--with-mpiexec' in help
configure=list(prior['configure_options'])
replace={'--prefix=':str(prefix),'PETSC_ARCH=':arch,'--with-cc=':str(mpi/'bin/mpicc'),'--with-cxx=':str(mpi/'bin/mpicxx'),'--with-blaslapack-lib=':'['+str(lib/'liblapack.so')+','+str(lib/'libblas.so')+']'}
for i,arg in enumerate(configure):
    for key,value in replace.items():
        if arg.startswith(key):configure[i]=key+value
configure.append('--with-mpiexec='+wrapper)
manifest['configure_command']=configure;write('cuda_build',manifest)
step(configure,'petsc_configure')
step(['make','-j8','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'petsc_make')
step(['make','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'install'],'petsc_install')
step(['make','PETSC_DIR='+str(prefix),'PETSC_ARCH=','check'],'petsc_check')
assert '#define PETSC_HAVE_CUDA 1' in (prefix/'include/petscconf.h').read_text()
linked=subprocess.check_output(['ldd',str(prefix/'lib/libpetsc.so')],env=ENV,text=True)
assert 'not found' not in linked
(L/'petsc_linked_libraries.log').write_text(linked)
manifest.update(status='PASS',self_test='PASS',library_sha256=hashlib.sha256((prefix/'lib/libpetsc.so').read_bytes()).hexdigest(),linked_libraries=linked)
write('cuda_build',manifest)

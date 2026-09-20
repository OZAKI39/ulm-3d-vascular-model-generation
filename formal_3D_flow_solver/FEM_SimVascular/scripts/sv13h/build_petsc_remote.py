"""Unpatched PETSc 3.19.6 with verified CUDA12 and preserved project MPI."""
import hashlib,json,shutil,tarfile
from pathlib import Path
from runner_remote import *
assert load('cuda12_runtime')['status']==load('mpi_with_cuda12')['status']=='PASS'
assert not (R/'petsc_cuda12_build.json').exists()
prior=json.loads((BASE/'configs/prior_cuda_build.json').read_text())
mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text())
cuda=Path(load('cuda12_install')['prefix']);mpiprefix=Path(mpi['prefix'])
source=BASE/'external/petsc-3.19.6-cuda12-clean';prefix=BASE/'external/petsc-cuda12';arch='arch-sv13h-cuda12'
archive=BASE.parent/'sv1_3/petsc-3.19.6.tar.gz'
assert digest(archive)==prior['source_archive_sha256'] and not source.exists()
original={}
with tarfile.open(archive) as t:
    for member in t:
        if member.isfile():original[str(Path(member.name).relative_to('petsc-3.19.6'))]=hashlib.sha256(t.extractfile(member).read()).hexdigest()
    t.extractall(BASE/'external',filter='data')
(BASE/'external/petsc-3.19.6').rename(source)
write('petsc_original_source_files',{'archive_sha256':digest(archive),'files':original})
libs=BASE/'external/system-lib';libs.mkdir()
for name in ('blas','lapack'):(libs/('lib'+name+'.so')).symlink_to('/usr/lib/x86_64-linux-gnu/lib'+name+'.so.3')
cuda_libraries=[cuda/'lib64'/('lib'+name+'.so') for name in ('cudart','cufft','cublas','cusparse','cusolver','curand','nvToolsExt')]
assert all(p.exists() for p in cuda_libraries),'CUDA12_LIBRARY_MISSING'
configure=['/usr/bin/python3','-B','configure','--prefix='+str(prefix),'PETSC_ARCH='+arch,'--with-cc='+str(mpiprefix/'bin/mpicc'),'--with-cxx='+str(mpiprefix/'bin/mpicxx'),'--with-mpiexec='+mpi['working_launcher'],'--with-fc=0','--with-fortran-bindings=0','--with-debugging=0','--with-shared-libraries=1','--with-scalar-type=real','--with-precision=double','--with-64-bit-indices=0','--with-blaslapack-lib=['+str(libs/'liblapack.so')+','+str(libs/'libblas.so')+']','--with-x=0','--with-cuda=1','--with-cudac='+str(cuda/'bin/nvcc'),'--with-cuda-include='+str(cuda/'include'),'--with-cuda-lib=['+','.join(map(str,cuda_libraries))+']','--with-cuda-arch=89','--with-cxx-dialect=C++17','--with-hip=0','--with-opencl=0']
manifest={'status':'BUILDING','version':prior['version'],'source_archive_sha256':digest(archive),'CUDA_release':load('cuda12_runtime')['selected_release'],'CUDA_prefix':str(cuda),'MPI_prefix':str(mpiprefix),'PETSC_DIR':str(source),'PETSC_ARCH':arch,'prefix':str(prefix),'configure_command':configure,'steps':[],'configure_exit':None,'make_exit':None,'source_unmodified':None,'cuda_enabled':None,'library_sha256':None}
write('petsc_cuda12_build',manifest)
def step(args,name):
    r=run(args,name,cwd=source,timeout=3600,extra_env={'PETSC_OPTIONS':'-use_gpu_aware_mpi 0'})
    manifest['steps'].append(r)
    if name=='petsc_cuda12_configure':
        manifest['configure_exit']=r['exit_code']
        if (source/'configure.log').exists():shutil.copyfile(source/'configure.log',L/'petsc_cuda12_configure_detail.log')
    if name=='petsc_cuda12_make':manifest['make_exit']=r['exit_code']
    write('petsc_cuda12_build',manifest)
    if not okay(r):
        manifest.update(status='FAIL',reason='PETSC_CUDA12_BUILD_FAIL',failed_step=name)
        if name=='petsc_cuda12_make' and any(x in text(r) for x in ("has no member named 'clockRate'","has no member named 'memoryClockRate'",'namespace "thrust" has no member class "unary_function"')):manifest['reason']='UNEXPECTED_CUDA12_API_FAILURE'
        write('petsc_cuda12_build',manifest);print(text(r)[-7000:],flush=True);raise SystemExit(1)
    return r
step(configure,'petsc_cuda12_configure')
conf=(source/arch/'include/petscconf.h').read_text();assert '#define PETSC_HAVE_CUDA 1' in conf
manifest['cuda_enabled']=True
step(['make','-j8','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'petsc_cuda12_make')
shutil.copyfile(source/arch/'lib/petsc/conf/make.log',L/'petsc_cuda12_make_detail.log')
step(['make','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'install'],'petsc_cuda12_install')
step(['make','PETSC_DIR='+str(prefix),'PETSC_ARCH=','check'],'petsc_cuda12_check')
changed=[name for name,h in original.items() if not (source/name).exists() or digest(source/name)!=h]
assert not changed,'PETSC_SOURCE_CHANGED'
ldd=run(['ldd',prefix/'lib/libpetsc.so'],'petsc_cuda12_ldd');assert okay(ldd) and 'not found' not in text(ldd)
compiler=run(['/usr/bin/g++','--version'],'petsc_cuda12_host_compiler')
manifest.update(status='PASS',source_unmodified=True,library_sha256=digest(prefix/'lib/libpetsc.so'),linked_libraries=text(ldd),compiler=text(compiler),self_test='PASS',build_time_s=sum(r['wall_time_s'] for r in manifest['steps']),source_files_verified=len(original))
write('petsc_cuda12_build',manifest)
print('PETSc CUDA12 configure/make/install/check PASS',flush=True)

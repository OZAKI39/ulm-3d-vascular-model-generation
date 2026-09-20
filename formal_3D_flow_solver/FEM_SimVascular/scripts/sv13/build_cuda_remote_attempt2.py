#!/usr/bin/env python3
"""Native PETSc 3.19.6 CUDA build in an isolated remote stage directory."""
import hashlib,json,os,shutil,subprocess,tarfile,time
from pathlib import Path
ROOT=Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3')
LOG=ROOT/'logs_attempt2';LOG.mkdir(exist_ok=True)
source=ROOT/'petsc-3.19.6';prefix=ROOT/'install/petsc-cuda'
archive=ROOT/'petsc-3.19.6.tar.gz'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest={'status':'BUILDING','version':'3.19.6','cuda_toolkit':'13.2','source_archive_sha256':digest(archive),
          'source':str(source),'prefix':str(prefix),'steps':[],'scientific_source_patched':False}
def record(): (ROOT/'cuda_build_manifest_attempt2.json').write_text(json.dumps(manifest,indent=2)+'\n')
assert not (ROOT/'cuda_build_manifest_attempt2.json').exists(),'Do not overwrite a build attempt'
if not source.exists():
    with tarfile.open(archive) as t:t.extractall(ROOT,filter='data')
lib=ROOT/'system-lib';lib.mkdir(exist_ok=True)
for name in ('blas','lapack'):
    p=lib/f'lib{name}.so'
    if not p.exists():p.symlink_to(f'/usr/lib/x86_64-linux-gnu/lib{name}.so.3')
env=dict(os.environ,PATH='/usr/bin:/bin:/usr/local/cuda/bin',OMP_NUM_THREADS='1',
         OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',CUDA_VISIBLE_DEVICES='0')
for k in ('PETSC_DIR','PETSC_ARCH','PETSC_OPTIONS','CC','CXX','FC','CFLAGS','CXXFLAGS','LDFLAGS','LD_LIBRARY_PATH'):env.pop(k,None)
def run(args,name):
    log=LOG/(name+'.log');start=time.monotonic()
    print('START '+name,flush=True)
    with log.open('x') as out:p=subprocess.run(args,cwd=source,env=env,stdout=out,stderr=subprocess.STDOUT)
    manifest['steps'].append({'command':args,'cwd':str(source),'exit_code':p.returncode,'elapsed_s':time.monotonic()-start,'log':str(log),'sha256':digest(log)})
    record();print('END '+name+' exit='+str(p.returncode),flush=True)
    if p.returncode:
        manifest.update(status='BLOCKED',reason='CUDA_PETSC_BUILD',failed_step=name);record()
        print(log.read_text(errors='replace')[-12000:],flush=True)
        raise SystemExit(p.returncode)
run(['/usr/bin/python3','-B','configure','--help'],'configure_help')
help_text=(LOG/'configure_help.log').read_text()
for option in ('--with-cuda','--with-cuda-dir','--with-cuda-arch','--with-cxx-dialect'):assert option in help_text
configure=['/usr/bin/python3','-B','configure','--prefix='+str(prefix),'PETSC_ARCH=arch-native-cuda-explicit',
    '--with-cc=/usr/bin/mpicc','--with-cxx=/usr/bin/mpicxx','--with-fc=0','--with-fortran-bindings=0',
    '--with-debugging=0','--with-shared-libraries=1','--with-scalar-type=real','--with-precision=double',
    '--with-64-bit-indices=0','--with-blaslapack-lib=['+str(lib/'liblapack.so')+','+str(lib/'libblas.so')+']',
    '--with-x=0','--with-cuda=1','--with-cudac=/usr/local/cuda/bin/nvcc','--with-cuda-include=/usr/local/cuda/include','--with-cuda-lib=['+','.join('/usr/local/cuda/lib64/lib'+name+'.so' for name in ('cudart','cufft','cublas','cusparse','cusolver','curand'))+']','--with-cuda-arch=89',
    '--with-cxx-dialect=C++17','--with-hip=0','--with-opencl=0']
manifest['configure_options']=configure;manifest['build_setup_correction']='Use documented explicit CUDA library list; missing legacy optional nvToolsExt library excluded; no source patch';record()
run(configure,'configure')
run(['/usr/bin/make','-j4','PETSC_DIR='+str(source),'PETSC_ARCH=arch-native-cuda-explicit','all'],'make')
run(['/usr/bin/make','PETSC_DIR='+str(source),'PETSC_ARCH=arch-native-cuda-explicit','install'],'install')
conf=(prefix/'include/petscconf.h').read_text();assert '#define PETSC_HAVE_CUDA 1' in conf
linked=subprocess.check_output(['ldd',str(prefix/'lib/libpetsc.so')],env=env,text=True)
(LOG/'linked_libraries.txt').write_text(linked)
manifest.update(status='PASS',cuda_support=True,binary_sha256=digest(prefix/'lib/libpetsc.so'),
                linked_libraries=linked,configuration_header_sha256=digest(prefix/'include/petscconf.h'))
record();print('CUDA PETSC BUILD PASS',flush=True)

#!/usr/bin/env python3
"""Build a minimal official PETSc and the unchanged pinned svMultiPhysics."""
import json, os, subprocess, sys, tarfile, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,write_json
LOG=ROOT/'logs/sv1_1'; REPORT=ROOT/'reports/sv1_1'
vtk=ROOT/'external/install/vtk'; prefix=ROOT/'external/install/sv1_1/petsc'
env=dict(os.environ,PATH='/usr/bin:/bin',LD_LIBRARY_PATH=f'{prefix}/lib:{vtk}/lib',
         PKG_CONFIG_PATH='',OMP_NUM_THREADS='1',CMAKE_BUILD_PARALLEL_LEVEL='4')
for key in ('PETSC_DIR','PETSC_ARCH','PETSC_OPTIONS','CC','CXX','FC','CFLAGS','CXXFLAGS','LDFLAGS'):
    env.pop(key,None)
manifest={'status':'BUILDING','steps':[],'version':'3.19.6','petsc_path':str(prefix)}
def run(cmd,name,cwd=ROOT):
    cmd=list(map(str,cmd)); start=time.monotonic(); log=LOG/(name+'.log'); resource=LOG/(name+'_resources.txt')
    if log.exists():
        suffix=1
        while log.with_suffix('.attempt'+str(suffix)+'.log').exists():suffix+=1
        log.rename(log.with_suffix('.attempt'+str(suffix)+'.log'))
        if resource.exists():resource.rename(resource.with_suffix('.attempt'+str(suffix)+'.txt'))
    print(name,flush=True)
    with log.open('w') as out:
        result=subprocess.run(['/usr/bin/time','-v','-o',str(resource),*cmd],cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT)
    meta={'command':cmd,'cwd':str(cwd),'exit_code':result.returncode,'elapsed_s':time.monotonic()-start,'log':str(log.relative_to(ROOT)),'resources':str(resource.relative_to(ROOT))}
    manifest['steps'].append(meta);write_json(REPORT/'petsc_build_manifest.json',manifest)
    if result.returncode:
        manifest.update(status='FAIL',reason='PETSC_BUILD_FAIL');write_json(REPORT/'petsc_build_manifest.json',manifest)
        print(log.read_text()[-8000:]);raise SystemExit(result.returncode)
    return meta

source=ROOT/'external/dependencies/sv1_1/petsc-3.19.6'
if not source.exists():
    source.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(ROOT/'external/downloads/sv1_1/petsc-3.19.6.tar.gz') as tar:tar.extractall(source.parent,filter='data')
system_lib=ROOT/'external/dependencies/sv1_1/system-lib';system_lib.mkdir(parents=True,exist_ok=True)
for name in ('blas','lapack'):
    link=system_lib/('lib'+name+'.so')
    if not link.exists():link.symlink_to('/lib/x86_64-linux-gnu/lib'+name+'.so.3')
if (source/'configure.log').exists():
    archive=LOG/('petsc_configure_detail_'+str(time.time_ns())+'.log')
    archive.write_bytes((source/'configure.log').read_bytes())
configure=['/usr/bin/python3','-B','configure','--prefix='+str(prefix),'PETSC_ARCH=arch-native-opt',
           '--with-cc=/usr/bin/mpicc','--with-cxx=/usr/bin/mpicxx','--with-fc=0',
           '--with-fortran-bindings=0','--with-debugging=0','--with-shared-libraries=1',
           '--with-scalar-type=real','--with-precision=double','--with-64-bit-indices=0',
           '--with-blaslapack-lib=['+str(system_lib/'liblapack.so')+','+str(system_lib/'libblas.so')+']',
           '--with-x=0','--with-cuda=0','--with-hip=0','--with-opencl=0']
run(configure,'petsc_configure',source)
assert (source/'arch-native-opt/include/petscconf.h').is_file(), 'PETSc configure did not generate configuration'
run(['/usr/bin/make','-j4','PETSC_DIR='+str(source),'PETSC_ARCH=arch-native-opt','all'],'petsc_build',source)
run(['/usr/bin/make','PETSC_DIR='+str(source),'PETSC_ARCH=arch-native-opt','install'],'petsc_install',source)
run(['/usr/bin/make','PETSC_DIR='+str(prefix),'PETSC_ARCH=','check'],'petsc_check',source)
solver=ROOT/'external/svMultiPhysics'; commit=subprocess.check_output(['git','-C',str(solver),'rev-parse','HEAD'],text=True).strip()
assert commit=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
assert not subprocess.check_output(['git','-C',str(solver),'status','--porcelain'],text=True).strip()
sv=ROOT/'external/SimVascularDistribution/usr/local/sv/simvascular/2023-05-31'
build=ROOT/'external/build/sv1_1/native-petsc'
extra=['-DSV_USE_TETGEN:BOOL=OFF','-DTETGEN_LIBRARY_NAME:FILEPATH='+str(sv/'lib/lib_simvascular_thirdparty_tetgen.a'),
       '-DCMAKE_CXX_FLAGS:STRING=-I'+str(sv/'include/thirdparty/tetgen'),
       '-DBLAS_LIBRARIES:FILEPATH=/lib/x86_64-linux-gnu/libblas.so.3','-DLAPACK_LIBRARIES:FILEPATH=/lib/x86_64-linux-gnu/liblapack.so.3']
run(['/usr/bin/cmake','-S',solver,'-B',build,'-DCMAKE_C_COMPILER=/usr/bin/gcc','-DCMAKE_CXX_COMPILER=/usr/bin/g++',
     '-DCMAKE_BUILD_TYPE=Release','-DSV_LOCAL_VTK_PATH='+str(vtk),'-DSV_PETSC_DIR='+str(prefix),
     '-DSV_ADDITIONAL_CMAKE_ARGS='+';'.join(extra)],'solver_petsc_configure')
run(['/usr/bin/cmake','--build',build,'--parallel','4'],'solver_petsc_build')
binary=build/'svMultiPhysics-build/bin/svmultiphysics'
linked=subprocess.check_output(['/usr/bin/ldd',str(binary)],env=env,text=True)
assert 'libpetsc' in linked and 'not found' not in linked
manifest.update(status='PASS',commit=commit,source_clean=True,executable=str(binary),executable_sha256=sha256(binary),
                linked_libraries=linked,linked_libraries_complete=True,runtime_library_path=env['LD_LIBRARY_PATH'],
                compiler=subprocess.check_output(['/usr/bin/g++','--version'],text=True),
                mpi=subprocess.check_output(['/usr/bin/mpiexec','--version'],text=True),
                petsc_header_sha256=sha256(prefix/'include/petscversion.h'),
                petsc_library_sha256=sha256(prefix/'lib/libpetsc.so'),
                mesh_library_recompiled=False,additional_algebra_packages=[])
write_json(REPORT/'petsc_build_manifest.json',manifest)
print('PETSC BUILD PASS',flush=True)

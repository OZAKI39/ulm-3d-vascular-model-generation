"""Build the unchanged pinned solver against the validated winner only."""
import json, re, shutil, subprocess, tarfile
from pathlib import Path
from runner_remote import *
assert load('petsc_gpu_smoke')['status']=='PASS'
assert load('ghost_probe_repair_03')['status']=='PASS' and load('ghost_coherence')['status']=='PASS'
w=load('compatibility_winner');prefix=Path(w['prefix']);wrapper=w['candidate_wrapper']
mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text());mpiprefix=Path(mpi['prefix'])
manifest=json.loads((BASE/'configs/solver_build_input_manifest.json').read_text())
archive=BASE.parent/'sv1_3j/external/solver_build_inputs.tar.gz';assert digest(archive)==manifest['sha256']
inputs=BASE/'external/solver-inputs';assert not inputs.exists();inputs.mkdir()
with tarfile.open(archive) as t:t.extractall(inputs,filter='data')
source=inputs/'svMultiPhysics';vtk=inputs/'vtk';tetgen=inputs/'tetgen'
head=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
assert head==manifest['source_git']['head']=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
names=subprocess.check_output(['git','-C',str(source),'ls-files','-z']).decode().split('\0')
original={n:digest(source/n) for n in names if n and (source/n).is_file()}
write('svmp_source_before',{'commit':head,'files':original,'source':str(source)})
build=BASE/'external/build/svmp_gpu_ghostfix';assert not build.exists()
extra_env={'OMPI_CC':'/usr/bin/gcc-12','OMPI_CXX':'/usr/bin/g++-12',
    'LD_LIBRARY_PATH':str(prefix/'lib')+':'+str(vtk/'lib'),'CMAKE_BUILD_PARALLEL_LEVEL':'4'}
args=['-DSV_USE_TETGEN:BOOL=OFF','-DTETGEN_LIBRARY_NAME:FILEPATH='+str(tetgen/'lib/lib_simvascular_thirdparty_tetgen.a'),
    '-DCMAKE_CXX_FLAGS:STRING=-I'+str(tetgen/'include'),'-DBLAS_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/libblas.so.3',
    '-DLAPACK_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/liblapack.so.3',
    '-DMPI_C_COMPILER:FILEPATH='+str(mpiprefix/'bin/mpicc'),'-DMPI_CXX_COMPILER:FILEPATH='+str(mpiprefix/'bin/mpicxx')]
cmd=['cmake','-S',source,'-B',build,'-DCMAKE_C_COMPILER=/usr/bin/gcc-12','-DCMAKE_CXX_COMPILER=/usr/bin/g++-12',
    '-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_STANDARD=17','-DSV_LOCAL_VTK_PATH='+str(vtk),
    '-DSV_PETSC_DIR='+str(prefix),'-DSV_ADDITIONAL_CMAKE_ARGS='+';'.join(args)]
d={'status':'BUILDING','source':str(source),'commit':head,'build':str(build),'PETSc_prefix':str(prefix),
    'MPI_prefix':str(mpiprefix),'VTK_prefix':str(vtk),'steps':[],'runtime_library_path':extra_env['LD_LIBRARY_PATH']}
def execute(args,name,timeout=3600):
    r=run(args,name,timeout=timeout,cuda=wrapper,extra_env=extra_env)
    d['steps'].append(r);write('svmp_gpu_build',d)
    if not okay(r):d.update(status='FAIL',reason='SVMULTIPHYSICS_GPU_BUILD_FAIL');write('svmp_gpu_build',d);raise SystemExit(1)
    return r
execute(cmd,'svmp_configure')
execute(['cmake','--build',build,'--parallel','4'],'svmp_make')
binary=build/'svMultiPhysics-build/bin/svmultiphysics';assert binary.is_file()
ldd=execute(['ldd',binary],'svmp_ldd');linked=text(ldd)
readelf=execute(['readelf','-d',binary],'svmp_readelf')
assert 'not found' not in linked
resolved=re.search(r'libpetsc[^\s]*\s+=>\s+(\S+)',linked)
assert resolved and Path(resolved[1]).resolve()==(prefix/'lib/libpetsc.so').resolve(),'SVMULTIPHYSICS_LINKED_WRONG_PETSC'
mpi_resolved=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',linked)
assert mpi_resolved and Path(mpi_resolved[1]).resolve()==(mpiprefix/'lib/libmpi.so').resolve(),'SVMULTIPHYSICS_LINKED_WRONG_MPI'
assert 'sv1_3g/external/gpu_mpi' not in linked
changed=[n for n,h in original.items() if not (source/n).is_file() or digest(source/n)!=h]
assert not changed,'SVMULTIPHYSICS_SOURCE_CHANGED'
after=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip();assert not after
evidence=R/'svmp_build_evidence';evidence.mkdir()
for p in (build/'CMakeCache.txt',build/'svMultiPhysics-build/CMakeCache.txt'):
    shutil.copyfile(p,evidence/(p.parent.name+'_CMakeCache.txt'))
d.update(status='PASS',executable=str(binary),executable_sha256=digest(binary),source_unmodified=True,
    source_files_verified=len(original),source_modifications=changed,source_git_status=after,linked_libraries=linked,
    linked_PETSc=resolved[1],PETSc_library_sha256=digest(prefix/'lib/libpetsc.so'))
write('svmp_gpu_build',d)
write('svmp_gpu_link',{'status':'PASS','resolved_PETSc':str(Path(resolved[1]).resolve()),'expected_PETSc':str((prefix/'lib/libpetsc.so').resolve()),
    'ldd':ldd,'readelf':readelf,'source_commit':head,'source_unmodified':True,'MPI_prefix':str(mpiprefix),'resolved_MPI':str(Path(mpi_resolved[1]).resolve()),'expected_MPI':str((mpiprefix/'lib/libmpi.so').resolve())})
print('Same-commit CUDA svMultiPhysics build and linkage PASS',flush=True)

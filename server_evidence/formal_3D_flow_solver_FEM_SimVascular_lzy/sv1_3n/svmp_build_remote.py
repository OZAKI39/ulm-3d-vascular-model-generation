"""Clean-build the pinned solver plus the recorded compatibility adapter."""
import json, re, shutil, subprocess, tarfile, sys
from pathlib import Path
from runner_remote import *
backend=sys.argv[1];attempt=sys.argv[2];assert backend in ('cpu','gpu13','gpu123','oldcompat')
if backend=='oldcompat':
 old=json.loads((BASE/'configs/baseline_L_compatibility_winner.json').read_text())
 old.update(version='3.19.6',reference_only=True)
 write('petsc_oldcompat_build',old)
w=load('petsc_'+backend+'_build');assert w['status']=='PASS'
prefix=Path(w['prefix']);wrapper=w['candidate_wrapper'];mpiprefix=Path(w['MPI_prefix'])
if backend not in ('cpu','oldcompat'):assert load('petsc_gpu_smoke')['status']=='PASS'
manifest=json.loads((BASE/'configs/baseline_L_solver_build_input_manifest.json').read_text())
archive=BASE.parent/'sv1_3j/external/solver_build_inputs.tar.gz';assert digest(archive)==manifest['sha256']
inputs=BASE/'external/solver-inputs'
if not inputs.exists():
 inputs.mkdir()
 with tarfile.open(archive) as t:t.extractall(inputs,filter='data')
baseline=inputs/'svMultiPhysics';source=BASE/'external'/('svmp_'+backend+'_'+attempt);assert not source.exists();shutil.copytree(baseline,source);vtk=inputs/'vtk';tetgen=inputs/'tetgen'
head=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
assert head==manifest['source_git']['head']=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
names=subprocess.check_output(['git','-C',str(source),'ls-files','-z']).decode().split('\0')
original={n:digest(source/n) for n in names if n and (source/n).is_file()}
write('svmp_'+backend+'_'+attempt+'_source_before',{'commit':head,'files':original,'source':str(source)})
patch=BASE/'patches/svmp_petsc325_compat.patch'
if attempt!='initial':
 assert patch.exists()
 applied=run(['patch','--batch','--fuzz=0','-p1','-i',patch],'svmp_'+backend+'_'+attempt+'_patch',cwd=source);assert okay(applied)
expected={n:digest(source/n) for n in original}
modified=[n for n in original if original[n]!=expected[n]]
write('svmp_'+backend+'_'+attempt+'_source_expected',dict(files=expected,modified=modified))
build=BASE/'external/build'/('svmp_petsc319_compat_cpu' if backend=='oldcompat' else 'svmp_petsc325_'+('cpu' if backend=='cpu' else 'gpu'))/attempt;assert not build.exists()
cc=w.get('CC','/usr/bin/gcc-12');cxx=w.get('CXX','/usr/bin/g++-12')
extra_env={'OMPI_CC':cc,'OMPI_CXX':cxx,
    'LD_LIBRARY_PATH':str(prefix/'lib')+':'+str(vtk/'lib'),'CMAKE_BUILD_PARALLEL_LEVEL':'4'}
args=['-DSV_USE_TETGEN:BOOL=OFF','-DTETGEN_LIBRARY_NAME:FILEPATH='+str(tetgen/'lib/lib_simvascular_thirdparty_tetgen.a'),
    '-DCMAKE_CXX_FLAGS:STRING=-I'+str(tetgen/'include'),'-DBLAS_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/libblas.so.3',
    '-DLAPACK_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/liblapack.so.3',
    '-DMPI_C_COMPILER:FILEPATH='+str(mpiprefix/'bin/mpicc'),'-DMPI_CXX_COMPILER:FILEPATH='+str(mpiprefix/'bin/mpicxx')]
cmd=['cmake','-S',source,'-B',build,'-DCMAKE_C_COMPILER='+cc,'-DCMAKE_CXX_COMPILER='+cxx,
    '-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_STANDARD=17','-DSV_LOCAL_VTK_PATH='+str(vtk),
    '-DSV_PETSC_DIR='+str(prefix),'-DSV_ADDITIONAL_CMAKE_ARGS='+';'.join(args)]
d={'status':'BUILDING','source':str(source),'commit':head,'build':str(build),'PETSc_prefix':str(prefix),
    'MPI_prefix':str(mpiprefix),'VTK_prefix':str(vtk),'steps':[],'runtime_library_path':extra_env['LD_LIBRARY_PATH']}
def execute(args,name,timeout=3600):
    r=run(args,'svmp_'+backend+'_'+attempt+'_'+name,timeout=timeout,cuda=wrapper,extra_env=extra_env)
    d['steps'].append(r);write('svmp_'+backend+'_build',d);write('svmp_'+backend+'_'+attempt+'_build',d)
    if not okay(r):d.update(status='FAIL',reason='SVMULTIPHYSICS_GPU_BUILD_FAIL');write('svmp_'+backend+'_build',d);write('svmp_'+backend+'_'+attempt+'_build',d);raise SystemExit(1)
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
changed=[n for n,h in expected.items() if not (source/n).is_file() or digest(source/n)!=h]
assert not changed,'SVMULTIPHYSICS_SOURCE_CHANGED'
after=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
evidence=R/('svmp_'+backend+'_'+attempt+'_build_evidence');evidence.mkdir()
for p in (build/'CMakeCache.txt',build/'svMultiPhysics-build/CMakeCache.txt'):
    shutil.copyfile(p,evidence/(p.parent.name+'_CMakeCache.txt'))
d.update(status='PASS',executable=str(binary),executable_sha256=digest(binary),source_unmodified=not modified,
    source_files_verified=len(original),source_modifications=modified,source_git_status=after,linked_libraries=linked,
    linked_PETSc=resolved[1],PETSc_library_sha256=digest(prefix/'lib/libpetsc.so'))
write('svmp_'+backend+'_build',d);write('svmp_'+backend+'_'+attempt+'_build',d)
write('svmp_'+backend+'_link',{'status':'PASS','resolved_PETSc':str(Path(resolved[1]).resolve()),'expected_PETSc':str((prefix/'lib/libpetsc.so').resolve()),
    'ldd':ldd,'readelf':readelf,'version':w['version'],'source_commit':head,'source_unmodified':not modified,'patch_sha256':digest(patch) if modified else None,'MPI_prefix':str(mpiprefix),'resolved_MPI':str(Path(mpi_resolved[1]).resolve()),'expected_MPI':str((mpiprefix/'lib/libmpi.so').resolve())})
print('Same-commit CUDA svMultiPhysics build and linkage PASS',flush=True)

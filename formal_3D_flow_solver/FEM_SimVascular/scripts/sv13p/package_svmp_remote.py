"""Build unchanged Stage O scientific source against an isolated GPU PETSc."""
import sys,re,shutil
from runner_remote import *
kind=sys.argv[1];assert kind in ('hypre','amgx')
ref=json.loads((BASE/'configs/reference_freeze.json').read_text());sv=ref['svmp'];w=load('petsc_'+kind+'_build');assert w['status']=='PASS'
patch=json.loads((BASE/'configs/source_patch.json').read_text())
source=BASE/('external/svMultiPhysics-'+kind);shutil.copytree(sv['source'],source)
for n,h in patch['before'].items():assert digest(source/n)==h,n
runtime=sv['runtime_library_path'].replace(ref['PETSc_build']['prefix'],w['prefix'])
build=BASE/('external/build/svmp_gpu_'+kind);cmd=sv['steps'][0]['command'][1:]
cmd=[str(source) if a==sv['source'] else str(build) if a==sv['build'] else a.replace(ref['PETSc_build']['prefix'],w['prefix']) for a in cmd]
extra={'LD_LIBRARY_PATH':runtime,'OMPI_CC':w['CC'],'OMPI_CXX':w['CXX'],'CMAKE_BUILD_PARALLEL_LEVEL':'4'}
steps=[]
for command,name in [(cmd,'configure'),(['cmake','--build',build,'--parallel','4'],'make')]:
 d=run(command,kind+'_svmp_'+name,timeout=3600,cuda=w['candidate_wrapper'],extra_env=extra);steps.append(d);assert okay(d)
binary=build/'svMultiPhysics-build/bin/svmultiphysics';assert binary.exists()
d=run(['ldd',binary],kind+'_svmp_ldd',cuda=w['candidate_wrapper'],extra_env=extra);steps.append(d);linked=text(d);assert okay(d) and 'not found' not in linked
p=re.search(r'libpetsc[^\s]*\s+=>\s+(\S+)',linked);m=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',linked)
assert p and Path(p[1]).resolve()==(Path(w['prefix'])/'lib/libpetsc.so').resolve()
assert m and Path(m[1]).resolve()==(Path(w['MPI_prefix'])/'lib/libmpi.so').resolve()
write('svmp_'+kind+'_build',dict(status='PASS',source=str(source),build=str(build),executable=str(binary),executable_sha256=digest(binary),runtime_library_path=runtime,PETSc_library_sha256=digest(p[1]),resolved_PETSc=str(Path(p[1]).resolve()),resolved_MPI=str(Path(m[1]).resolve()),PETSc_build=w,steps=steps,source_files_verified=len(patch['before'])))
print(kind+' isolated solver build PASS.',flush=True)

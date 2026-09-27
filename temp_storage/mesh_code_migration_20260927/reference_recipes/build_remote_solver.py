"""Snapshot protected stacks; build only the isolated P1 adapter."""
import re,shutil
from runner_remote import *
from environment_remote import snapshot
ref=json.loads((BASE/'configs/reference_freeze.json').read_text());patch=json.loads((BASE/'configs/source_patch.json').read_text());sv=ref['svmp'];w=ref['PETSc_build']
assert digest(sv['executable'])==ref['GPU_baseline']['solver_sha256']
assert digest(Path(w['prefix'])/'lib/libpetsc.so')==ref['GPU_baseline']['PETSc_library_sha256']
assert not (R/'pre_install_environment.json').exists();write('pre_install_environment',snapshot())
O=BASE.parent/'sv1_3p';out=BASE/'outputs';out.mkdir(exist_ok=True)
shutil.copytree(BASE.parent/'sv1_3o/outputs/SV_MESH',out/'SV_MESH');shutil.copyfile(O/'outputs/REAL_VASCULAR_GPU_PC_WINNER/solver.xml',out/'vascular_template.xml')
source=BASE/'external/svMultiPhysics-reuse';shutil.copytree(sv['source'],source)
for n,h in patch['before'].items():assert digest(source/n)==h,n
p=run(['patch','--batch','--fuzz=0','-p1','-i',BASE/'patches/ilu_rebuild_policy.patch'],'reuse_patch',cwd=source);assert okay(p)
for n,h in patch['after'].items():assert digest(source/n)==h,n
build=BASE/'external/build/svmp_gpu_reuse';cmd=sv['steps'][0]['command'][1:]
cmd=[str(source) if a==sv['source'] else str(build) if a==sv['build'] else a for a in cmd]
extra={'LD_LIBRARY_PATH':sv['runtime_library_path'],'OMPI_CC':w['CC'],'OMPI_CXX':w['CXX'],'CMAKE_BUILD_PARALLEL_LEVEL':'4'}
steps=[]
for command,name in [(cmd,'configure'),(['cmake','--build',build,'--parallel','4'],'make')]:
 d=run(command,'reuse_'+name,timeout=3600,cuda=w['candidate_wrapper'],extra_env=extra);steps.append(d);assert okay(d)
binary=build/'svMultiPhysics-build/bin/svmultiphysics';assert binary.exists()
d=run(['ldd',binary],'reuse_ldd',cuda=w['candidate_wrapper'],extra_env=extra);steps.append(d);linked=text(d);assert okay(d) and 'not found' not in linked
p=re.search(r'libpetsc[^\s]*\s+=>\s+(\S+)',linked);m=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',linked)
assert p and Path(p[1]).resolve()==(Path(w['prefix'])/'lib/libpetsc.so').resolve()
assert m and Path(m[1]).resolve()==(Path(w['MPI_prefix'])/'lib/libmpi.so').resolve()
write('svmp_reuse_build',dict(status='PASS',source=str(source),build=str(build),executable=str(binary),executable_sha256=digest(binary),runtime_library_path=sv['runtime_library_path'],PETSc_library_sha256=digest(p[1]),resolved_PETSc=str(Path(p[1]).resolve()),resolved_MPI=str(Path(m[1]).resolve()),PETSc_build=w,steps=steps,source_files_verified=len(patch['after'])))
print('Stage Q isolated solver build PASS.',flush=True)

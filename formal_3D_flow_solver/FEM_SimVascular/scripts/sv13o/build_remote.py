"""Read-only baseline audit, stage-local scientific input copies and one clean build."""
import json,re,shutil,subprocess
from pathlib import Path
from runner_remote import *
from environment_remote import snapshot
ref=json.loads((BASE/'configs/reference_freeze.json').read_text())
patch=json.loads((BASE/'configs/source_patch.json').read_text())
policy=json.loads((BASE/'configs/policy.json').read_text())
N=BASE.parent/'sv1_3n';sv=ref['svmp'];w=ref['PETSc_build']
assert digest(sv['executable'])==ref['GPU_baseline']['solver_sha256']
assert digest(Path(w['prefix'])/'lib/libpetsc.so')==ref['GPU_baseline']['PETSc_library_sha256']
assert digest(N/'outputs/REAL_VASCULAR_GPU/1-procs/stFile_060.bin')==policy['checkpoint_sha256']
assert (R/'pre_install_environment.json').exists()
source=BASE/'external/svMultiPhysics'
import tarfile
with tarfile.open(BASE/'missing_tracked_source.tar.gz') as t:
 for member in t.getmembers():assert not (source/member.name).exists(),member.name
 t.extractall(source,filter='data')
for n,h in patch['before'].items():assert digest(source/n)==h,n
assert digest(BASE/'patches/final_output_on_stop.patch')==patch['patch_sha256']
p=run(['patch','--batch','--fuzz=0','-p1','-i',BASE/'patches/final_output_on_stop.patch'],'output_patch',cwd=source);assert okay(p)
for n,h in patch['after'].items():assert digest(source/n)==h,n
build=BASE/'external/build';assert not build.exists()
cmd=sv['steps'][0]['command'][1:]
cmd=[str(source) if a==sv['source'] else str(build) if a==sv['build'] else a for a in cmd]
extra={'LD_LIBRARY_PATH':sv['runtime_library_path'],'OMPI_CC':w['CC'],'OMPI_CXX':w['CXX'],'CMAKE_BUILD_PARALLEL_LEVEL':'4'}
steps=[]
for command,name in [(cmd,'configure'),(['cmake','--build',build,'--parallel','4'],'make')]:
 d=run(command,'svmp_'+name,timeout=3600,cuda=w['candidate_wrapper'],extra_env=extra);steps.append(d);assert okay(d)
binary=build/'svMultiPhysics-build/bin/svmultiphysics';assert binary.exists()
ldd=run(['ldd',binary],'svmp_ldd',cuda=w['candidate_wrapper'],extra_env=extra);steps.append(ldd)
linked=text(ldd);assert okay(ldd) and 'not found' not in linked
p=re.search(r'libpetsc[^\s]*\s+=>\s+(\S+)',linked);m=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',linked)
assert p and Path(p[1]).resolve()==(Path(w['prefix'])/'lib/libpetsc.so').resolve()
assert m and Path(m[1]).resolve()==(Path(w['MPI_prefix'])/'lib/libmpi.so').resolve()
for n,h in patch['after'].items():assert digest(source/n)==h,n
write('svmp_build',dict(status='PASS',source=str(source),build=str(build),executable=str(binary),executable_sha256=digest(binary),runtime_library_path=sv['runtime_library_path'],PETSc_library_sha256=digest(p[1]),resolved_PETSc=str(Path(p[1]).resolve()),resolved_MPI=str(Path(m[1]).resolve()),steps=steps,source_files_verified=len(patch['after']),only_Stage_O_source_change=patch['changed_files'],Stage_N_lifecycle_preserved=True))
print('Stage O clean GPU build PASS; frozen PETSc/MPI link and output-only source change verified.',flush=True)

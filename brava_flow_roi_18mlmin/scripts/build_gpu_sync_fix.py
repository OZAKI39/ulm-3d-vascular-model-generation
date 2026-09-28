"""Isolated one-object rebuild; original solver, build and headers are never edited."""
from pathlib import Path
import json,shutil,hashlib,subprocess,shlex,difflib
ROOT=Path('/workspace/brava_flow_roi_18mlmin_20260928');out=ROOT/'gpu_solver_fix';out.mkdir(exist_ok=True)
build=json.loads(Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports/svmp_reuse_build.json').read_text())
source=Path(build['source']);b=Path(build['build'])/'svMultiPhysics-build/Source/solver';objdir=b/'CMakeFiles/svmultiphysics.dir'
snapshot=out/'source'
if not snapshot.exists():shutil.copytree(source,snapshot,ignore=shutil.ignore_patterns('.git'))
assert (snapshot/'Code/Source/solver/petsc_impl.cpp').read_bytes()==(source/'Code/Source/solver/petsc_impl.cpp').read_bytes()
p=snapshot/'Code/Source/solver/petsc_impl.cpp';before=p.read_text()
needle='    /* Fill the ghost vertices with correct values. */'
patch='''    // The ghost-local sequential view shares the host buffer, but its state
    // does not track writes made to the parent CUDA vector by KSPSolve.
    // Synchronize owned solution entries before the CPU ghost-local export.
    const PetscScalar* brava_owned_host = nullptr;
    PetscCallAbort(MPI_COMM_WORLD, VecGetArrayRead(psol[cEq].b, &brava_owned_host));
    PetscCallAbort(MPI_COMM_WORLD, VecRestoreArrayRead(psol[cEq].b, &brava_owned_host));

'''
assert before.count(needle)>=1;p.write_text(before.replace(needle,patch+needle,1))
(out/'gpu_host_sync.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),p.read_text().splitlines(True),fromfile='original/petsc_impl.cpp',tofile='patched/petsc_impl.cpp')))
flags={}
for line in (objdir/'flags.make').read_text().splitlines():
 if ' = ' in line:k,v=line.split(' = ',1);flags[k]=shlex.split(v)
obj=out/'PetscLinearAlgebra.cpp.o';exe=out/'svmultiphysics'
compile=['/usr/bin/g++',*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],'-c',str(snapshot/'Code/Source/solver/PetscLinearAlgebra.cpp'),'-o',str(obj)]
link=shlex.split((objdir/'link.txt').read_text());idx=link.index('CMakeFiles/svmultiphysics.dir/PetscLinearAlgebra.cpp.o');link[idx]=str(obj);link[link.index('-o')+1]=str(exe)
with (out/'build.log').open('w') as f:
 subprocess.run(compile,cwd=b,stdout=f,stderr=subprocess.STDOUT,check=True)
 subprocess.run(link,cwd=b,stdout=f,stderr=subprocess.STDOUT,check=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=dict(status='BUILT_NOT_YET_VALIDATED',base_executable=build['executable'],base_executable_sha256=sha(build['executable']),executable=str(exe),executable_sha256=sha(exe),patch_sha256=sha(out/'gpu_host_sync.patch'),
 source_modified=['Code/Source/solver/petsc_impl.cpp'],original_source_sha256=sha(source/'Code/Source/solver/petsc_impl.cpp'),patched_source_sha256=sha(p),
 scope='Explicit device-to-host solution synchronization only; no FEM, mesh, boundary, constitutive, time integration or residual tolerance change',compile_command=compile,link_command=link,link_cwd=str(b))
(out/'build.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:v for k,v in manifest.items() if k not in ['compile_command','link_command']},indent=2))

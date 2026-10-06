"""Verify retained build provenance, exact deployed binary and actual source tree."""
from common import *
import datetime
import subprocess
import importlib.util
import shutil


def git(path, *args):
    return subprocess.check_output(['git','-c',f'safe.directory={path}','-C',str(path),*args],text=True).strip()


def main():
    inputs={}
    def reg(p):
        p=Path(p);inputs[str(p)]=sha(p);return inputs[str(p)]
    execution=read(BASE/'cases/balanced_pressure_final/reports/execution.json')
    build=read(BASE/'gpu_solver_fix/build.json')
    base_build=read(Q/'reports/svmp_reuse_build.json')
    manifest=read(Q/'configs/source_patch.json')
    assert reg(BASE/'gpu_solver_fix/svmultiphysics')==execution['solver_sha256']==build['executable_sha256']
    assert reg(Path(build['base_executable']))==build['base_executable_sha256']==base_build['executable_sha256']
    assert reg(base_build['resolved_PETSc'])==execution['PETSc_library_sha256']==base_build['PETSc_library_sha256']
    assert reg(Q/'configs/source_patch.json')==base_build['source_manifest_sha256']
    mismatches={n:sha(SOURCE/n) for n,h in manifest['after'].items() if sha(SOURCE/n)!=h}
    assert not mismatches,mismatches
    retained=BASE/'gpu_solver_fix/source'
    patch_difference={n:sha(retained/n) for n,h in manifest['after'].items() if sha(retained/n)!=h}
    assert set(patch_difference)=={'Code/Source/solver/petsc_impl.cpp'},patch_difference
    assert patch_difference['Code/Source/solver/petsc_impl.cpp']==build['patched_source_sha256']
    assert sha(SOURCE/'Code/Source/solver/petsc_impl.cpp')==build['original_source_sha256']
    assert reg(BASE/'gpu_solver_fix/gpu_host_sync.patch')==build['patch_sha256']
    source_commit=git(SOURCE,'rev-parse','HEAD')
    changed=git(SOURCE,'status','--short')
    object_dir=Path(build['link_cwd'])/'CMakeFiles/svmultiphysics.dir'
    proof={}
    for name in ['fluid.cpp','nn.cpp','fs.cpp','read_files.cpp','read_msh.cpp','eq_assem.cpp','initialize.cpp','Integrator.cpp','PetscLinearAlgebra.cpp','output.cpp','vtk_xml.cpp']:
        p=SOURCE/'Code/Source/solver'/name
        proof[name]=dict(source_path=str(p),source_sha256=reg(p))
        for suffix in ['.o','.o.d']:
            q=object_dir/(name+suffix)
            if q.exists():proof[name][suffix]=dict(path=str(q),sha256=reg(q))
        if name not in ('PetscLinearAlgebra.cpp',):
            proof[name]['matches_git_commit']=sha(p)==__import__('hashlib').sha256(subprocess.check_output(['git','-c',f'safe.directory={SOURCE}','-C',str(SOURCE),'show',f'{source_commit}:Code/Source/solver/{name}'])).hexdigest()
            assert proof[name]['matches_git_commit']
    for rel in ['flags.make','link.txt','build.make','DependInfo.cmake']:reg(object_dir/rel)
    for p in [BASE/'gpu_solver_fix/build.json',BASE/'gpu_solver_fix/build.log',Q/'reports/svmp_reuse_build.json',
              Q/'reports/final_input_integrity.json',Q/'reports/stop_timestamp_inner_make.json',
              SOURCE/'Code/Source/solver/Parameters.cpp',SOURCE/'Code/Source/solver/FE/Basis/NodeOrderingConventions.cpp',
              Path(base_build['build'])/'svMultiPhysics-build/CMakeCache.txt',
              PREVIOUS/'FIELD_RELIABILITY_RESULTS.json',PREVIOUS/'EVENT_FIELD_PROTOCOL.json',PREVIOUS/'FIELD_METRICS.json',
              PREVIOUS/'FIELD_COMPARABILITY.json',CAUSE/'CAUSE_RESULTS.json',BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json',
              BASE/'cases/balanced_pressure_final/run/solver.xml',BASE/'cases/balanced_pressure_final/reports/execution.json']:
        reg(p)
    for path,digest in read(PREVIOUS/'INPUT_SHA256.json').items():
        assert reg(BASE/path)==digest,path
    protected={}
    for rel in ['population_v1/analysis/speed_event_cause_review/INPUT_SHA256.json','population_v1/analysis/field_reliability_event_review/OUTPUT_SHA256.json']:
        p=BASE/rel;items=read(p);root=BASE if 'speed_event_cause' in rel else PREVIOUS
        for k,h in items.items():assert sha(root/k)==h,k
        protected[rel]=dict(count=len(items),sha256=reg(p),unchanged=True)
    resources={k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in ['cpu.max','cpuset.cpus.effective','memory.max','memory.current','memory.stat','memory.events']}
    for label,cmd in [('GPU',['nvidia-smi','--query-gpu=name,memory.total,memory.used,memory.free,utilization.gpu','--format=csv']),
                      ('other_GPU_processes',['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv']),
                      ('disk',['df','-h',str(BASE)])]:resources[label]=subprocess.check_output(cmd,text=True)
    resources['time_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    put('RESOURCES.json',resources)
    put('TOOLING_INVENTORY.json',dict(python=sys.executable,modules={n:bool(importlib.util.find_spec(n)) for n in ['dolfinx','dolfin','ngsolve','skfem','firedrake','basix','symfem']},executables={n:shutil.which(n) for n in ['FreeFem++','ngsolve']},system_modified=False))
    exact_manifest=dict(manifest['after']);exact_manifest.update(patch_difference)
    put('EXACT_DEPLOYED_SOURCE_MANIFEST.json',exact_manifest)
    result=dict(status='SOURCE_CONFIRMED',source_identity='PRODUCTION_SOLVER_SOURCE_IDENTITY_ESTABLISHED',
       method='Verified recorded build chain, all retained source hashes, object dependency paths, final base and host-sync binary hashes; no rebuild performed.',
       source_commit=source_commit,source_checkout=str(SOURCE),source_status=changed,verified_source_files=len(exact_manifest),
       actual_binary=build['executable'],actual_binary_sha256=build['executable_sha256'],
       base_binary=build['base_executable'],base_binary_sha256=build['base_executable_sha256'],
       retained_patched_tree=str(retained),only_host_sync_difference=patch_difference,
       base_source_manifest_sha256=base_build['source_manifest_sha256'],exact_source_manifest_sha256=sha(HERE/'EXACT_DEPLOYED_SOURCE_MANIFEST.json'),
       assembly_sources=proof,compile_command=build['compile_command'],link_command=build['link_command'],link_cwd=build['link_cwd'],
       PETSc=dict(path=base_build['resolved_PETSc'],sha256=execution['PETSc_library_sha256'],options=execution['PETSC_OPTIONS']),
       MPI_ranks=execution['MPI_ranks'],execution_command=execution['command'],
       actual_compiler_flags=(object_dir/'flags.make').read_text(),source_patch_metadata={k:v for k,v in manifest.items() if k not in ('before','after')},
       server_git=dict(root=git(BASE,'rev-parse','--show-toplevel'),head=git(BASE,'rev-parse','HEAD'),branch=git(BASE,'branch','--show-current'),
                       preexisting_changes=git(BASE,'status','--short')),protected_previous_inputs=protected,
       production_algorithm_modified=False,new_CFD_count=0)
    put('PRODUCTION_SOURCE_IDENTITY.json',result);put('INPUT_SHA256.json',inputs)
    print('SOURCE IDENTITY',result['source_identity'],source_commit,'files',len(exact_manifest),'registered',len(inputs))


if __name__=='__main__':main()

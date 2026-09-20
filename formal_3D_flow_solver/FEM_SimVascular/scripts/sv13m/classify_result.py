"""Record the measured Stage M outcome without promoting partial flow completion."""
import json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3m';sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13m import ghost_gate
from sv_validation.provenance import sha256
load=lambda n:json.loads((R/(n+'.json')).read_text())
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
for n in ('compatibility_winner','petsc_gpu_smoke','petsc_cuda_types','svmp_gpu_build','svmp_gpu_link','patch_integrity','remote_preservation','native_artifacts'):
 source=R/'remote'/(n+'.json')
 if source.exists():shutil.copyfile(source,R/(n+'.json'))
shutil.copyfile(R/'remote/repair_03_petsc_build.json',R/'petsc_build.json')
ghost_gate(load('ghost_probe_after'));assert load('ghost_coherence')['status']=='PASS'
iterations=[]
audit=load('upstream_audit')
r1={k:audit[k] for k in ('modified_files','modified_functions','added_lines','deleted_lines')}
r1.update(iteration='repair_01',trigger='Baseline CUDA ghost rejects derived types',root_cause='Concrete equality misses seqcuda/mpicuda',upstream_commit=audit['first_fixed_commit'],first_fixed_tag=audit['first_fixed_tag'],tests_added=['Four CPU/CUDA ghost variants, each three repetitions'],before_result='CPU PASS / CUDA seq and MPI FAIL',after_result='CUDA seq and CPU PASS; CUDA MPI still loses ghost metadata',scientific_impact='NONE — compatibility only');write('repair_01',r1)
for i in range(1,4):
 r=load(f'repair_0{i}')
 if i==2:r['after_result']='CPU and CUDA seq PASS; CUDA MPI forward PASS / reverse data FAIL'
 if i==3:r['after_result']='All four variants 3/3 PASS; all forward/reverse/local data correct; CPU/CUDA coherence and duplicate 3/3 PASS'
 r['patch_path']=f'patches/sv1_3m/repair_0{i}/petsc319_cuda_ghost_backport.patch';r['patch_sha256']=sha256(ROOT/r['patch_path'])
 write(f'repair_0{i}',r);iterations.append(r)
write('repair_iterations',dict(status='PASS',iterations=iterations,used=3,maximum=3,remaining=0))
smoke=load('svmp_gpu_smoke');assert smoke['status']=='FAIL' and smoke['MPI_error_detected']
reason='Official GPU smoke exit failure after MPI_Finalize; all three authorized compatibility repairs consumed. No additional source repair performed.'
for name in ('GPU_PROOF_20_acceptance','CPU_PROOF_20_1R_acceptance','science_equivalence','gpu_residency','ghost_transfer','benchmark'):
 write(name,dict(status='NOT_RUN',reason=reason,measurements=None,scientific_inputs_changed=False))
write('stage_result',dict(status='FAIL',classification='SVMULTIPHYSICS_GPU_SMOKE_FAIL',detail='POST_MPI_FINALIZE_GPU_RESOURCE_TEARDOWN_FAILURE',ghost_repair='PASS',standalone_GPU='PASS',svmp_build='PASS',official_smoke='FAIL',repair_iterations_used=3,repair_limit_reached=True,remaining_work_requires_new_authorization=True,vascular_proof='NOT_RUN',science_equivalence='NOT_RUN',GPU_residency='NOT_RUN',benchmark='NOT_RUN',production='CPU_EARLY_STOP_PRODUCTION',production_changed=False,full_GPU_steady_started=False,mesh_convergence_started=False,stack_label='PETSc 3.19.6 + upstream CUDA ghost compatibility backport',coherence_patch_note='Includes explicitly labeled local compatibility adaptation, not all hunks are verbatim upstream.'))
E=R/'exit_failure_evidence';E.mkdir(exist_ok=True)
for file in ('main.cpp','PetscLinearAlgebra.cpp','petsc_impl.cpp'):
 p=ROOT/'external/svMultiPhysics/Code/Source/solver'/file;shutil.copyfile(p,E/file)
write('exit_failure_diagnosis',dict(status='DIAGNOSED',official_exit_code=smoke['exit_code'],original_ghost_error_absent=True,MPI_ERR_TYPE_absent=True,new_error='MPI_Comm_rank() function was called after MPI_FINALIZE was invoked',official_linear_solves=smoke['linear_solves'],official_nonlinear_steps=smoke['steps_completed'],native_VTU_count=smoke['VTU_count'],native_fields_finite=smoke['velocity_finite'] and smoke['pressure_finite'],fresh_reload=smoke['reload_pass'],source_findings=['PetscLinearAlgebra::finalize() is empty (line 93)','main.cpp calls that no-op then MPI_Finalize (lines 676-683)','petsc_destroy_all contains PetscFinalize but is not called on this normal exit path'],gdb_findings='Separate diagnostic replay catches SIGSEGV in PETSc CUDA event_pool global destructor / PoolAllocator during __cxa_finalize; not a KSP or ghost update frame.',scope_limit='3/3 repairs used; solver source must remain unchanged. No fourth repair or error suppression.',evidence=['logs/sv1_3m/remote/official_fluid_gpu_smoke.log','logs/sv1_3m/official_exit_diagnostic.log','reports/sv1_3m/exit_failure_evidence/'],root_cause_confidence='Observed exit lifecycle and GPU resource teardown failure; exact allocator corruption causality not proven by a memory sanitizer.',official_contract='https://petsc.org/release/manualpages/Sys/PetscFinalize/'))
print(json.dumps(load('stage_result'),indent=2))

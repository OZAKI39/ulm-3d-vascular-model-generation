from runner_remote import *
sv=load('svmp_reuse_build');w=sv['PETSc_build'];prefix=Path(w['prefix']);header=Path(sv['source'])/'Code/Source/solver/sv13q_reuse.h'
extra={'LD_LIBRARY_PATH':sv['runtime_library_path'],'OMPI_CXX':w['CXX']}
binary=BASE/'benchmarks/recovery_probe'
compile=run([Path(w['MPI_prefix'])/'bin/mpicxx','-std=c++17','-O2','-I'+str(header.parent),'-I'+str(prefix/'include'),BASE/'recovery_probe.cpp','-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lpetsc','-o',binary],'recovery_probe_compile_v3',cuda=w['candidate_wrapper'],extra_env=extra);assert okay(compile)
N=BASE.parent/'sv1_3n';launcher=json.loads((N/'configs/baseline_L_mpi_application_gate.json').read_text())['working_launcher']
r=run([launcher,'-n','1',w['candidate_wrapper'],binary,'-skip_petscrc','-use_gpu_aware_mpi','0','-ksp_converged_reason','-log_view',':recovery_probe_profile.txt','-log_view_gpu_time'],'recovery_probe_run_v3',cuda=w['candidate_wrapper'],extra_env=extra);assert okay(r)
s=text(r);assert 'recovery=STALE_ILU_RECOVERED' in s and 'no_third_attempt=1' in s
write('adaptive_recovery_probe',dict(status='PASS',synthetic_fault_injection=True,not_CFD=True,not_performance_evidence=True,production_header_sha256=digest(header),compile=compile,run=r,original_RHS_restored=True,same_operator=True,zero_initial_guess=True,fresh_failure_retained=True,at_most_one_retry=True,initial_probe_issues='First standalone run omitted frozen MPI option and stopped before solving. Second passed recovery assertions but harness cleanup used old PETSc context-destruction API. Both failed logs retained. Only standalone invocation/cleanup corrected; production header unchanged.'))

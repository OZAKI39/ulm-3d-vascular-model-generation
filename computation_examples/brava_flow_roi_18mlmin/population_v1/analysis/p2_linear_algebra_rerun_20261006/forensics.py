from common import *
import re,shutil,datetime
r=read(HERE/'P2_LINEAR_SYSTEM_RESOURCE_RESULT.json');p=HERE/'monitoring/p2_system/job.log';log=p.read_text()
pre=read(HERE/'P2_SYSTEM_PRELAUNCH.json');samples=r['resource_samples'];late=samples[-5:]
pressure_nnz=[int(x) for x in re.findall(r'rows=259330, cols=259330\n\s*total: nonzeros=(\d+)',log)]
put('P2_PRESSURE_MATRIX_OBSERVATIONS.json',dict(source=str(p.relative_to(HERE)),log_sha256=sha(p),observed_pressure_matrix_NNZ=pressure_nnz,maximum_observed_NNZ=max(pressure_nnz),matches_boolean_graph_bound=max(pressure_nnz)==read(HERE/'CPU_SCHUR_GRAPH_BUDGET.json')['two_hop_node_nnz'],warning='Repeated MatView summaries do not identify distinct simultaneously live allocations; do not sum them as memory copies.'))
put('P2_STOP_DIAGNOSIS.json',dict(status='SUPERVISED_RESOURCE_STOP_NO_OOM',exact_trigger=r['stop_reason'],predeclared_tree_RSS_limit_GiB=24,sampled_tree_RSS_peak_GiB=r['actual_P2_job_peak_tree_RSS_GiB'],host_forecast_GiB=pre['host_job_forecast_GiB'],observed_peak_minus_forecast_GiB=r['actual_P2_job_peak_tree_RSS_GiB']-pre['host_job_forecast_GiB'],headroom_before_stop_GiB=r['minimum_conservative_headroom_GiB'],required_headroom_GiB=r['required_headroom_GiB'],instantaneous_headroom_threshold_failed=False,post_setup_headroom_test='NOT_EVALUABLE_SETUP_INCOMPLETE',late_growth_GiB_per_s=r['last20sample_growth_GiB_per_s'],last_five_samples=late,exact_allocating_function='NOT_CAPTURED; do not infer an exact call stack from matrix summaries',observed_stage='Inside first PC setup: velocity hierarchy summaries reach 6x6; then pressure 259330x259330 sparse product/submatrix summaries; no iteration-0 after_PC_setup callback',pressure_product_graph_bound=read(HERE/'CPU_SCHUR_GRAPH_BUDGET.json'),precision_of_peak='Sampled process-tree peak at interrupted setup, LOWER BOUND on a completed setup requirement; per-rank records are earlier high-water lower bounds',CUDA_error=False,host_OOM_delta=r['host_oom_event_delta'],no_parameter_retry=True,forecast_correction='Previous 17.93GiB CPU envelope underestimated actual setup; it is invalid for future admission. Existing P2_SYSTEM_ATTEMPT guard refuses rerun.',remaining_unknowns=['Exact CPU PC transient/factor allocation attribution','Completed PC setup peak','P2 linear convergence','Steady P2 physical field']))
# Keep concise source-backed evidence for the mature backend, not an installed library copy.
old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');src=Path(old['PETSc']['path']).parent.parent.parent/'source_gpu13_configfix';spec={'src/ksp/pc/impls/fieldsplit/fieldsplit.c':[(685,776),(853,918)],'src/ksp/ksp/utils/schurm/schurm.c':[(782,855)],'src/mat/impls/aij/mpi/mpimatmatmult.c':[(799,851),(901,941)],'src/ksp/pc/impls/gamg/gamg.c':[(130,152)]};out=[];meta={}
for n,ranges in spec.items():
 f=src/n;lines=f.read_text().splitlines();meta[n]=dict(path=str(f),sha256=sha(f),line_ranges=ranges);out.append(n)
 for a,b in ranges:out.extend(f'{i+1}: {lines[i]}' for i in range(a-1,min(b,len(lines))))
(HERE/'CPU_PC_NATIVE_SOURCE_EVIDENCE.txt').write_text('\n'.join(out)+'\n');put('CPU_PC_SOURCE_IDENTITY.json',meta)
put('ATTEMPT_HISTORY.json',dict(preparation_error='JSON inheritance assumed dict; numeric-count list wrapped. No scientific operation involved.',P1_attempt1='Monitoring assertion inspected implicit Schur operator instead of its CPU preconditioning matrix. Stopped at first monitor, no accepted time step.',observer_fix='KSPGetOperators(...,nullptr,&P) for AMG levels. Candidate options/science/tolerances/RA policy unchanged.',P1_attempt2='PASS, 9 nonlinear iterations and one completed restart segment step10->11',P2_attempt1='RSS guard stopped setup; zero P2 linear iterations, zero CFD steps, no retry',allocation_smoke='Deliberately oversized request denied BEFORE CUDA allocation; expected exit86, not OOM'))
# Preserve full receipts/logs, small native observers, and guard logs as evidence.
evidence=HERE/'evidence';evidence.mkdir(exist_ok=True)
receipts=[]
for d in sorted((HERE/'monitoring').iterdir()):
 rp=d/'EXECUTION.json'
 if not rp.exists():continue
 rec=read(rp);receipts.append(dict(name=d.name,exit_code=rec['exit_code'],stop_reason=rec['stop_reason'],elapsed_s=rec['elapsed_s'],peak_tree_RSS_GiB=rec['peak_tree_RSS_GiB'],receipt=f'evidence/{d.name}_EXECUTION.json'))
 shutil.copy2(rp,evidence/f'{d.name}_EXECUTION.json');shutil.copy2(d/'job.log',evidence/f'{d.name}_job.log')
for case in ['p1_candidate','p1_candidate_fixed_observer','p2_system']:
 run=HERE/'runtime'/case/'run'
 for f in list(run.glob('algebra_rank_*.jsonl'))+list(run.glob('numeric_rank_*.json'))+[run/'cuda_admission.jsonl']:
  if f.exists():shutil.copy2(f,evidence/f'{case}_{f.name}')
put('EXECUTION_RECEIPT.json',dict(status='EXITED',jobs=receipts,total_supervised_job_seconds=sum(x['elapsed_s'] for x in receipts),total_is_not_development_wall_time=True,completed_P1_parity_steps=1,P2_single_system_attempts=1,P2_completed_linear_solves=0,P2_CFD_steps=0,other_project_jobs_modified=False,automatic_next_job=False,UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
put('RESULT_IDENTITY.json',dict(task='P2/P1 algebra rerun with CPU Schur PC and CUDA outer operator',authority='Current isolated rerun; prior result unchanged',candidate_plan_sha256=sha(HERE/'LINEAR_SOLVER_CANDIDATE_PLAN.json'),build=read(HERE/'BUILD_IDENTITY.json'),mesh_sha256=read(TARGET/'P2_CURVED_MESH_QA.json').get('curved_mesh_sha256','dbad7f2be2e4eb4a718f2db275b6ab3bf663915f229088367530db705407d48e'),inherited_manifest_sha256=sha(HERE/'INHERITED_EVIDENCE.json'),P1_baseline_restart=read(HERE/'P1_RESTART_CORRESPONDENCE.json'),production_modified=False))
print('Forensics and receipts saved',len(receipts))

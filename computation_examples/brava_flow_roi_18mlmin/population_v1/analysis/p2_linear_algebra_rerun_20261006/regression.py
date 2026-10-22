from common import *
import ast,subprocess,re,xml.etree.ElementTree as ET
from PIL import Image
checks=[]
def check(name,value):
 checks.append(dict(name=name,PASS=bool(value)))
p=read(HERE/'P1_LINEAR_SOLVER_PARITY.json');r=read(HERE/'P2_LINEAR_SYSTEM_RESOURCE_RESULT.json');d=read(HERE/'P2_DOF_STORAGE_AUDIT.json');g=read(HERE/'P2_LINEAR_ALGEBRA_DECISION_GATE.json');plan=read(HERE/'LINEAR_SOLVER_CANDIDATE_PLAN.json')
check('one_predeclared_candidate',plan['candidate_count']==1)
check('P1_same_CSR_and_RHS',p['candidate']['matrix_identity']['CSR_bitwise_equal'] and p['candidate']['matrix_identity']['RHS_bitwise_equal'])
check('P1_matrix_roundoff_bound',p['candidate']['matrix_identity']['row_roundoff_units']<=128)
check('P1_solution_parity',p['candidate']['PASS'])
check('P1_CPU_hierarchy',p['candidate']['CPU_PC_hierarchy_verified'])
check('P1_native_CPU_copy_checks',p['candidate']['CPU_Pmat_identity_checks']==9)
check('P1_convergence_cost_predeclared',p['candidate']['total_KSP_iterations']<=10*p['baseline']['total_KSP_iterations'])
check('P1_XML_same_as_verified_segment',sha(HERE/'runtime/p1_candidate_fixed_observer/run/solver.xml')==sha(PRIOR/'runtime/p1_schur/run/solver.xml'))
check('P2_XML_same_science',sha(HERE/'runtime/p2_system/run/solver.xml')==sha(PRIOR/'runtime/p2_system/run/solver.xml'))
for n in ['inlet_exact_velocity.dat','outlet_exact_traction.vtp']:
 check('P2_BC_byte_identity_'+n,sha(HERE/'runtime/p2_system/run'/n)==sha(PRIOR/'runtime/p2_system/run'/n))
check('pressure_edge_count_identity',d['total_equation_rows']-d['mathematical_rows']==d['edge_pressure_storage_only_slots']==223119)
check('numeric_plus_zero_equals_structure',d['numeric_nnz']+d['numeric_structural_zero_entries']==d['scalar_structural_nnz']==111793312)
check('compact_reference_unchanged',sha(HERE/'COMPACT_MIXED_GRAPH_REFERENCE.npz')==sha(PRIOR/'COMPACT_MIXED_GRAPH_REFERENCE.npz'))
check('common_geometry_object_unchanged',sha(PREVIOUS/'runtime/patched.o')==read(HERE/'BUILD_IDENTITY.json')['fluid_object_sha256'])
check('P2_stopped_for_declared_RSS',r['stop_reason']=='PROCESS_TREE_RSS_LIMIT' and r['actual_P2_job_peak_tree_RSS_GiB']>24)
check('instant_headroom_not_misreported',r['sampled_headroom_at_all_recorded_points_ge20percent'] and not r['headroom_20percent_proven'])
check('continued_growth_recorded',not r['non_growing_gate'] and r['last20sample_growth_GiB_per_s']>.005)
check('P2_setup_not_claimed_complete',r['setup']=='ATTEMPTED_NO_COMPLETION' and not r['CPU_PC_hierarchy_verified'])
check('P2_no_linear_result',r['single_system_result'] is None and not (HERE/'runtime/p2_system/run/single_result.json').exists())
check('no_cuda_allocation_failure',not r['CUDA_allocation_failure'] and r['CUDA_admission_denials']==0)
check('no_host_OOM',r['host_oom_event_delta']==0 and r['host_oom_kill_event_delta']==0)
check('guard_native_smoke',read(HERE/'CUDA_ADMISSION_GUARD.json')['status']=='PASS')
check('NO_is_bounded_scope',g['P2_PIPE_LINEAR_ALGEBRA_FEASIBLE']=='NO' and 'does not prove' in g['linear_feasibility_scope'])
check('high_order_output_still_FAIL',read(HERE/'HIGH_ORDER_OUTPUT_GATE.json')['status']=='FAIL' and not g['high_order_output_gate'])
check('full_CFD_not_allowed',not g['full_P2_CFD_allowed'] and g['full_P2_CFD_count']==0)
check('no_BraVa_particles',all(g[x]==0 for x in ['BraVa_CFD_count','trajectory_count','population_count','particle_dt_sweep_count']))
for name in ['P2_PIPE_RUN_IDENTITY.json','P2_PIPE_FIELD_VALIDATION.json','P2_PIPE_ERROR_DECOMPOSITION.json']:
 check(name+'_NOT_RUN',read(HERE/name)['status']=='NOT_RUN')
check('compact_not_implemented',read(HERE/'COMPACT_PROTOTYPE_GATE.json')['status']=='NOT_IMPLEMENTED')
for name,h in read(HERE/'INPUT_SHA256.json').items():check('protected:'+name,sha(BASE/name)==h)
for name,h in read(HERE/'INHERITED_EVIDENCE.json').items():check('prior_immutable:'+name,sha(BASE/name)==h)
# One-attempt gate test executes only the early refusal, no subprocess solver.
run=subprocess.run([sys.executable,str(HERE/'p2_system.py')],text=True,capture_output=True)
check('P2_rerun_refused_before_execution',run.returncode!=0 and 'P2_SYSTEM_ALREADY_ATTEMPTED_NO_RETRY' in run.stderr)
put('REENTRY_GUARD_TEST.json',dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr,scientific_operations=0))
for f in HERE.glob('*.py'):
 ast.parse(f.read_text());check('syntax:'+f.name,True)
for name in ['01_P2_MEMORY_BREAKDOWN.png','02_ASM_OVERLAP_AND_FILL.png','03_P2_SOLVER_RESOURCE_COMPARISON.png']:
 with Image.open(HERE/name) as im:check(name+'_300dpi',all(abs(x-300)<.1 for x in im.info['dpi']) and im.width>=3000)
for label,url in re.findall(r'\[([^\]]+)\]\(([^)]+)\)',(HERE/'P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md').read_text()):
 if url.startswith('http') or url in {'TEST_RESULTS.json','FINAL_INTEGRITY.json'}:continue
 check('report_link:'+label,(HERE/url).exists())
check('TASK_removed',not (HERE/'TASK.json').exists())
put('TEST_RESULTS.json',dict(status='PASS' if all(x['PASS'] for x in checks) else 'FAIL',checks=checks,count=len(checks),passed=sum(x['PASS'] for x in checks),scientific_precision_certification=False,scope='Integrity, regression and gate behavior; resource feasibility and output gates remain failed'))
print(len(checks),sum(x['PASS'] for x in checks));assert all(x['PASS'] for x in checks),[x for x in checks if not x['PASS']]

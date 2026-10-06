from common import *
p=read(HERE/'P1_LINEAR_SOLVER_PARITY.json');f=read(HERE/'P2_LINEAR_SYSTEM_RESOURCE_RESULT.json');o=read(HERE/'HIGH_ORDER_OUTPUT_GATE.json') if (HERE/'HIGH_ORDER_OUTPUT_GATE.json').exists() else {'status':'NOT_RUN'}
# Scientific field result is never inferred from a successful algebraic solve.
feasible='YES' if f['status']=='P2_LINEAR_SYSTEM_RESOURCE_PASS' else 'NO' if f['status']=='P2_LINEAR_SYSTEM_MEMORY_UNSAFE' else 'NOT_EVALUABLE'
if feasible=='YES':recommendation='P2P1_LINEAR_ALGEBRA_FEASIBLE_BUT_PIPE_NOT_RUN';reason='HIGH_ORDER_OUTPUT_GATE_FAILED' if o['status']=='FAIL' else 'FULL_CFD_PENDING'
elif feasible=='NO':recommendation='STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT';reason=f['status']
else:recommendation='P2P1_PIPE_RESULT_NOT_EVALUABLE';reason='NO_SAFE_LOW_MEMORY_LINEAR_SOLVER_CANDIDATE' if p['status']!='PASS' else f['status']
gates=dict(frozen_BC_geometry=read(TARGET/'P2_PIPE_READINESS_GATE.json')['nine_implementation_prerequisites_pass'],P1_parity=p['status']=='PASS',actual_P2_system_resource_margin=feasible=='YES',P2_linear_convergence=f.get('solve')=='CONVERGED',high_order_output=o['status']=='PASS')
put('P2_LINEAR_ALGEBRA_DECISION_GATE.json',dict(P2_PIPE_LINEAR_ALGEBRA_FEASIBLE=feasible,NEXT_RECOMMENDATION=recommendation,stop_reason=reason,gates=gates,full_CFD_authorized_by_gates=all(gates.values()),full_P2_CFD_runs=0,BraVa_P2_CFD=0,trajectory=0,population=0,particle_dt_sweeps=0,actual_P2_field_improvement='NOT_RUN',BraVa_pilot='NOT_EVALUABLE',compact_numbering_prototype='NOT_IMPLEMENTED',compact_reason='Only structural payload savings established; no demonstrated >=15% whole-job peak improvement. Native numbering unchanged.',next_resource_direction='HIGHER_HOST_MEMORY_SERVER' if feasible=='NO' else None))
for name in ['P2_PIPE_RUN_IDENTITY.json','P2_PIPE_FIELD_VALIDATION.json','P2_PIPE_ERROR_DECOMPOSITION.json']:
 put(name,dict(status='NOT_RUN',reason=reason,actual_P2_steady_field=False,CFD_runs=0,scientific_inference='NONE'))
print(feasible,recommendation,reason)

c=read(HERE/'COMPACT_MIXED_GRAPH_STATS.json');peak=f.get('actual_P2_job_peak_tree_RSS_GiB')
put('COMPACT_PROTOTYPE_GATE.json',dict(status='NOT_IMPLEMENTED',structural_root_cause_established=True,threshold_whole_job_peak_benefit=.15,known_one_CSR_payload_saving_GiB=c['CSR_value_index_saving_bytes']/2**30,known_CSR_saving_over_observed_P2_peak=c['CSR_value_index_saving_bytes']/2**30/peak if peak else None,whole_job_peak_saving_proven=False,PC_fill_or_AMG_change_proven=False,native_Val_would_remain_in_local_PETSc_only_prototype=True,eligibility='No demonstrated >=15% full-process peak improvement; mathematical compact graph is provided but no solver numbering change',P1_numbering_unchanged=True,P2_mathematical_DOF_unchanged=True))

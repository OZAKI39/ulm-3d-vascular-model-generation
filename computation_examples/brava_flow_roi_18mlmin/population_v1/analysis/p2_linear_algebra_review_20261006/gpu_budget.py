"""Fail closed before an allocation whose temporary workspace has no bound."""
def gate(budget, free_MiB, ranks=8, reserve_MiB=2048):
 if not budget or budget.get('status')!='BOUNDED' or budget.get('MPI_ranks')!=ranks:
  return dict(PASS=False,reason='GPU_MEMORY_PREFLIGHT_NOT_ESTABLISHED')
 if budget.get('source')!='native_workspace_query_without_numeric_allocation':
  return dict(PASS=False,reason='UNVERIFIED_WORKSPACE_BOUND')
 workspace=budget.get('sum_rank_simultaneous_workspace_bound_MiB');resident=budget.get('resident_increment_bound_MiB')
 if not isinstance(workspace,(int,float)) or not isinstance(resident,(int,float)) or min(workspace,resident)<0:
  return dict(PASS=False,reason='INVALID_WORKSPACE_BOUND')
 required=workspace+resident+reserve_MiB
 return dict(PASS=required<=free_MiB,reason='PASS' if required<=free_MiB else 'GPU_BUDGET_EXCEEDS_FREE',required_MiB=required,free_MiB=free_MiB,reserve_MiB=reserve_MiB)

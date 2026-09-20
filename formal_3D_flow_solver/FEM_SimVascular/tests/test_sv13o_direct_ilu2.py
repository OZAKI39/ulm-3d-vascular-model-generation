from sv13o_support import *
def test_actual_ASM_topology_precedes_B():
 t=accepted('asm_topology');assert set(t['blocks_per_solve'])=={1}
 assert t['local_rows']==t['global_rows']==281452
def test_B_changes_only_wrapper():
 d=result('PERF_B');assert d['status'] in ('PASS','FAIL')
 plan=json.loads((C/'runplans/PERF_B.json').read_text())
 assert '-pc_type ilu' in plan['PETSC_OPTIONS'] and '-pc_factor_levels 2' in plan['PETSC_OPTIONS']
 assert '-pc_type asm' not in plan['PETSC_OPTIONS']
 if d['status']=='PASS':
  for s in d['runtime_semantics']:assert s['PC']=='ilu' and s['sub_PC'] is None and s['fill_level']==2 and s['restart']==100
 else:assert d['errors']


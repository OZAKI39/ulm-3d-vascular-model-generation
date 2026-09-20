from sv13o_support import *
def test_C_is_conditional_and_preserves_other_settings():
 d=result('PERF_C');a=accepted('PERF_A_acceptance');b=result('PERF_B')
 if d['status']=='NOT_RUN':assert b['status']=='PASS' and 1-b['wall_time_s']/a['wall_time_s']>=.10;return
 assert b['status']!='PASS' or 1-b['wall_time_s']/a['wall_time_s']<.10
 plan=json.loads((C/'runplans/PERF_C.json').read_text())
 assert plan['PETSC_OPTIONS'].replace('-pc_factor_levels 1','-pc_factor_levels 2')==b['PETSC_OPTIONS']
 if d['status']=='PASS':assert all(s['fill_level']==1 and s['PC']=='ilu' for s in d['runtime_semantics'])
 else:assert d['status']=='FAIL' and d['errors']

def test_failed_C_log_still_proves_actual_ILU1_configuration():
 from sv_validation.sv13o import parse_runtime_semantics
 d=result('PERF_C');rows=parse_runtime_semantics((ROOT/'logs/sv1_3o/remote/PERF_C.log').read_text())
 assert rows and rows[0]['fill_level']==1
 for s in rows:semantics_gate(s,d['PETSC_OPTIONS'])
 assert d['status']=='FAIL' and d['KSP_reason_values']['DIVERGED_BREAKDOWN']==-5

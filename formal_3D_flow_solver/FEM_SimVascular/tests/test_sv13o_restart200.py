from sv13o_support import *
def test_D_uses_incumbent_PC_and_only_changes_restart():
 plan=json.loads((C/'candidate_D.json').read_text());base=result('PERF_'+plan['base_candidate']);d=result('PERF_D')
 assert plan['PETSC_OPTIONS']==base['PETSC_OPTIONS'].replace('-ksp_gmres_restart 100','-ksp_gmres_restart 200')
 if d['status']=='PASS':assert all(s['restart']==200 for s in d['runtime_semantics'])
 else:assert d['status']=='FAIL' and d['errors']
 assert all('-ksp_gmres_restart '+str(n) not in plan['PETSC_OPTIONS'] for n in (150,250,300,400))


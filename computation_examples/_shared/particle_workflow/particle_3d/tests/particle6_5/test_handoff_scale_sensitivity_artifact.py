import json,numpy as np

def test_all_floors_same_setup(report):
 d=json.loads((report/'data/08_sensitivity.json').read_text());rows=d['rows']
 assert len(rows)==9 and d['scientific_sensitivity_review']=='PENDING_USER_REVIEW'
 assert d['no_arbitrary_scientific_pass_threshold'] and d['all_finite'] and d['all_continuously_safe']
 assert {r['h_molecular_floor_m'] for r in rows}=={1.5e-9,2e-9,3e-9}
 assert len({r['horizon_s'] for r in rows})==1
 for div in [1,2,4]:
  subset=[r for r in rows if r['dt_divisor']==div];assert len({r['dt_s'] for r in subset})==1
  assert all(r['parameter_role']=='SENSITIVITY_ONLY' for r in subset)

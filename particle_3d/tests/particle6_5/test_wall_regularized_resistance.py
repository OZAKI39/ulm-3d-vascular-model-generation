import numpy as np

def test_wall_scalar_balance(scans):
 for r in scans[0]:
  expected=r['w']*6*np.pi*.00345312*r['radius_i_m']**2/max(r['h_geom_m'],2e-9)
  assert np.isclose(r['coefficient_kg_s'],expected,rtol=2e-15,atol=0)
  assert r['scalar_solution_error']<8*np.finfo(float).eps

def test_raw_gap_not_capped(scans):
 rows=[r for r in scans[0] if r['h_geom_m']<2e-9]
 assert rows and all(r['h_eff_m']==2e-9 and r['gap_m']==r['h_geom_m'] and r['interaction_state']=='BELOW_HANDOFF_TRIAL_REJECTED' for r in rows)

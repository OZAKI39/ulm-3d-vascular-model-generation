import numpy as np

def test_pairs_static(scans):
 for r in scans[1]:
  reff=r['radius_i_m']*r['radius_j_m']/(r['radius_i_m']+r['radius_j_m'])
  assert np.isclose(r['coefficient_kg_s'],r['w']*6*np.pi*.00345312*reff**2/r['h_eff_m'],rtol=2e-15,atol=0)
  assert r['scalar_solution_error']<8*np.finfo(float).eps
  if r['chi']>=.05:assert r['coefficient_kg_s']==0
  if r['h_geom_m']>r['h_lower_m'] and r['chi']<=.01:assert r['coefficient_kg_s']==r['old_p5_leading_uncapped_zeta_kg_s']

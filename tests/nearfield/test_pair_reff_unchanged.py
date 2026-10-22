def test_reff(scans):
 for r in scans[1]:
  a,b=r['radius_i_m'],r['radius_j_m'];assert r['resistance_radius_m']==a*b/(a+b)
  assert r['resistance_radius_m']!=r['a_ref_m']

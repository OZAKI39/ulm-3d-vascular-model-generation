def test_min_radius(scans):
 assert {r['radius_ratio'] for r in scans[1]}=={1,2,4}
 assert all(r['a_ref_m']==min(r['radius_i_m'],r['radius_j_m']) for r in scans[1])

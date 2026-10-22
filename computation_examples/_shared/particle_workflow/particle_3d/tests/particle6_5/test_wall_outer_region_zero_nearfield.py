import numpy as np

def test_outer(scans):
 rows=[r for r in scans[0] if r['chi']>=.05]
 assert rows and all(r['coefficient_kg_s']==0 and abs(r['v1_unconstrained_vn_ratio']-1)<=8*np.finfo(float).eps for r in rows)

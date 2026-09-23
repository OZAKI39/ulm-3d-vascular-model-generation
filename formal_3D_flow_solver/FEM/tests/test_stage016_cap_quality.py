from stage016_helpers import *
import pytest
from fem3d.cap_remesh import triangle_quality
from fem3d.planar_port import cap_gates

@pytest.mark.parametrize('name',CANDIDATES)
def test_actual_3d_cap_quality_and_interior_planarity(name):
 d=candidate(name);q,_=triangle_quality(d['points_m'],d['triangles'][d['facet_tags']!=1])
 assert np.isfinite(q).all() and np.count_nonzero(q<.1)==0 and np.quantile(q,.05)>=.45 and np.median(q)>=.70
 for r in audit_candidate(name).values():assert r['new_interior_max_plane_deviation_m']<=1e-15

def test_bad_quality_cannot_pass_a_sparse_budget():
 counts={'inlet':56,'outlet_01':44,'outlet_02':42,'outlet_03':49}
 for q in (np.full(191,.09),np.full(191,.4),np.full(191,.6)):
  assert cap_gates(counts,q,POLICY)['status']=='FAIL'

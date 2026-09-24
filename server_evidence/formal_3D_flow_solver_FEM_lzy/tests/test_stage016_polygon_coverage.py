from stage016_helpers import *
import pytest
from fem3d.cap_remesh import coverage_check

@pytest.mark.parametrize('name',CANDIDATES)
def test_actual_cap_is_one_nonoverlapping_polygon_disk(name):
 for r in audit_candidate(name).values():
  assert r['coverage']['holes']==r['coverage']['overlaps']==r['coverage']['triangles_outside_polygon']==0

def test_malicious_hole_rejected():
 p,t=square()
 with pytest.raises(ValueError):coverage_check(p[:,:2],t[1:],np.arange(4))

def test_malicious_overlapping_triangle_rejected():
 p,t=square()
 with pytest.raises(ValueError):coverage_check(p[:,:2],np.vstack([t,t[:1]]),np.arange(4))

def test_outside_polygon_even_if_rim_unchanged_rejected():
 p,t=square();p[4,:2]=[1.4,.5]
 with pytest.raises(ValueError):coverage_check(p[:,:2],t,np.arange(4))

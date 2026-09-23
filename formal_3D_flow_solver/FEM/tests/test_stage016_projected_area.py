from stage016_helpers import *
import pytest
from fem3d.planar_port import verify_metrics
from copy import deepcopy

@pytest.mark.parametrize('name',CANDIDATES)
def test_two_independent_projected_area_implementations(name):
 for result in audit_candidate(name).values():
  assert result['relative_projected_area_error']<=1e-12
  assert result['independent_projected_area_relative_error']<=1e-12

def test_malicious_changed_projected_area_rejected():
 p=CONTRACT['ports']['inlet'];m=deepcopy(p);m['formal_projected_area_m2']*=1+1e-8
 with pytest.raises(ValueError,match='Projected area'):verify_metrics(p,m,POLICY['geometry'])

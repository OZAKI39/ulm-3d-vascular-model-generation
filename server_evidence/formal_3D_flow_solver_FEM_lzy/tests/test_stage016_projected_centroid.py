from stage016_helpers import *
import pytest
from fem3d.planar_port import verify_metrics
from copy import deepcopy

@pytest.mark.parametrize('name',CANDIDATES)
def test_projected_polygon_and_triangle_centroids_agree(name):
 for r in audit_candidate(name).values():
  assert r['projected_centroid_displacement_m']<=1e-12
  assert r['independent_projected_centroid_error_m']<=1e-12

def test_shifted_centroid_rejected():
 p=CONTRACT['ports']['inlet'];m=deepcopy(p);m['projected_centroid_m'][0]+=1e-10
 with pytest.raises(ValueError,match='centroid'):verify_metrics(p,m,POLICY['geometry'])

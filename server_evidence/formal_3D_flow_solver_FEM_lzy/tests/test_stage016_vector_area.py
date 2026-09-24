from stage016_helpers import *
import pytest
from fem3d.planar_port import verify_metrics,polygon_metrics
from copy import deepcopy

@pytest.mark.parametrize('name',CANDIDATES)
def test_vector_is_rim_contour_invariant(name):
 for n,result in audit_candidate(name).items():
  p=CONTRACT['ports'][n]
  independent=np.cross(np.array(p['rim_coordinates_m'])-p['plane_origin_m'],np.roll(np.array(p['rim_coordinates_m'])-p['plane_origin_m'],-1,axis=0)).sum(axis=0)/2
  np.testing.assert_allclose(independent,result['A_vector_candidate'],rtol=1e-12,atol=1e-28)
  assert result['relative_vector_area_error']<=1e-12 and result['normal_dot_source']>=1-1e-12

def test_malicious_reversed_vector_rejected():
 p=CONTRACT['ports']['inlet'];m=deepcopy(p);m['vector_area_m2']=(-np.array(m['vector_area_m2'])).tolist()
 with pytest.raises(ValueError,match='Vector'):verify_metrics(p,m,POLICY['geometry'])

def test_reversed_rim_orientation_rejected():
 p=CONTRACT['ports']['inlet']
 with pytest.raises(ValueError,match='Reversed'):polygon_metrics(p['rim_coordinates_m'][::-1],p['plane_origin_m'],p['basis'])

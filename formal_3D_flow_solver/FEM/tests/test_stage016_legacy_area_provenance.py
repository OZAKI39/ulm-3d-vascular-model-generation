from stage016_helpers import *
import pytest
from fem3d.mesh_qc import triangle_geometry

@pytest.mark.parametrize('name',list(CONTRACT['ports']))
def test_legacy_area_and_centroid_preserved(name):
 p=CONTRACT['ports'][name];s=source();a,c,_=triangle_geometry(s['points_m'],s['triangles'][s['facet_tags']==p['entity_id']])
 assert float(a.sum())==p['legacy_scalar_area_m2']==p['legacy_scalar_surface_area_m2']
 np.testing.assert_array_equal(np.sum(c*a[:,None],axis=0)/a.sum(),p['legacy_scalar_area_centroid_m'])
 assert p['relative_difference_legacy_scalar_vs_projected']>1e-12

@pytest.mark.parametrize('candidate_name',CANDIDATES)
def test_scalar_difference_is_disclosed_but_does_not_reject_geometry(candidate_name):
 for r in audit_candidate(candidate_name).values():
  assert abs(r['candidate_scalar_vs_legacy_relative_difference'])>1e-12
  assert r['status']=='PASS' and r['legacy_scalar_area_is_gate'] is False

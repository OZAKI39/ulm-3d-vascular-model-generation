from stage017_helpers import *
import pytest
@pytest.mark.parametrize('shape',['ellipse','mild_irregular_convex_polygon'])
def test_actual_scaled_Gmsh_searches_match(shape):
 d=read(OUT/'synthetic/results.json');assert d['status']=='PASS';cases=d['shapes'][shape]['cases'];ref=cases[0]
 assert [c['scale'] for c in cases]==[.5,1.,4.]
 for c in cases[1:]:
  assert c['search']['selected_trial']==ref['search']['selected_trial'] and c['triangles']==ref['triangles']
  assert c['normalized_H']==pytest.approx(ref['normalized_H'])
  np.testing.assert_allclose(c['normalized_points'],ref['normalized_points'],atol=1e-10,rtol=0)
  for a,b in zip(c['search']['trials'],ref['search']['trials']):np.testing.assert_allclose(list(a['quality']['q_tri'].values()),list(b['quality']['q_tri'].values()),atol=1e-10,rtol=0)

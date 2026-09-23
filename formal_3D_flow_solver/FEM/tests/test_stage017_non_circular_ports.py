from stage017_helpers import *
import pytest
from fem3d.cap_remesh import coverage_check,cross2

@pytest.mark.parametrize('shape',['ellipse','mild_irregular_convex_polygon'])
def test_actual_non_circle_rim_is_preserved_and_covered(shape):
 record=read(OUT/'synthetic/results.json')['shapes'][shape];assert not record['circle_fit'] and record['actual_rim_used']
 for case in record['cases']:
  d=np.load(OUT/'synthetic'/shape/str(case['scale'])/(case['search']['selected_trial']+'.npz'));rim=d['xy'][:48]
  assert np.ptp(np.linalg.norm(rim,axis=1))>0
  assert coverage_check(d['xy'],d['triangles'],np.arange(48))['status']=='PASS'

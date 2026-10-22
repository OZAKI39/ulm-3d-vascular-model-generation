import numpy as np
import pytest
from network_1d0d.idealized_h0 import solve_operating_point,external_radius_variant

@pytest.mark.parametrize('factor',[.90,.95,1.,1.05,1.10])
def test_actual_global_scaled_radius_solve(domain,baseline,factor):
    result=solve_operating_point(domain,radii=domain.radius_m*factor)
    np.testing.assert_allclose(result['fractions'],baseline['fractions'],rtol=1e-10,atol=1e-12)
    assert abs(result['scaling_lambda']/baseline['scaling_lambda']-factor**-4)<1e-9

def test_exterior_perturbation_keeps_roi_geometry(domain):
    r,mask=external_radius_variant(domain,'O2',.95)
    assert mask.sum()==76
    inside=np.unique(domain.edges[domain.internal_roi])
    np.testing.assert_array_equal(r[inside],domain.radius_m[inside])
    with pytest.raises(ValueError,match='no saved downstream'):
        external_radius_variant(domain,'O3',.95)

import numpy as np
from particle_3d.validation_cases import AFFINE_ATOL, AFFINE_RTOL


def test_analytic_curl_and_symmetric_rate(affine):
    for name, expected in [("vorticity_s_inv", affine["vorticity"]), ("strain_rate_s_inv", affine["strain"])]:
        actual = getattr(affine["batch"], name)
        np.testing.assert_allclose(actual, np.broadcast_to(expected, actual.shape), atol=AFFINE_ATOL, rtol=AFFINE_RTOL)

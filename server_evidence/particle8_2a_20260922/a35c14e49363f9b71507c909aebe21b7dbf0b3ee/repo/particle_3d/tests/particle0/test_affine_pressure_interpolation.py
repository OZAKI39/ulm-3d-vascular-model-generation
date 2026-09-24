import numpy as np
from particle_3d.validation_cases import AFFINE_ATOL, AFFINE_RTOL


def test_affine_pressure_all_point_kinds(affine):
    print("max affine pressure error (Pa)", np.abs(affine["batch"].pressure_pa - affine["pressure"]).max())
    np.testing.assert_allclose(affine["batch"].pressure_pa, affine["pressure"], atol=AFFINE_ATOL, rtol=AFFINE_RTOL)

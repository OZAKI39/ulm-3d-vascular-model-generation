import numpy as np
from particle_3d.validation_cases import AFFINE_ATOL, AFFINE_RTOL


def test_affine_velocity_all_point_kinds(affine):
    error = np.abs(affine["batch"].velocity_m_s - affine["velocity"])
    print("max affine velocity error (m/s)", error.max(), "atol", AFFINE_ATOL, "rtol", AFFINE_RTOL)
    np.testing.assert_allclose(affine["batch"].velocity_m_s, affine["velocity"], atol=AFFINE_ATOL, rtol=AFFINE_RTOL)

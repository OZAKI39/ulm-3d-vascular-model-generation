import numpy as np
from particle_3d.validation_cases import AFFINE_ATOL, AFFINE_RTOL


def test_affine_gradient_index_convention(affine):
    expected = np.broadcast_to(affine["gradient"], affine["batch"].velocity_gradient_s_inv.shape)
    print("max affine gradient error (1/s)", np.abs(affine["batch"].velocity_gradient_s_inv - expected).max())
    np.testing.assert_allclose(affine["batch"].velocity_gradient_s_inv, expected, atol=AFFINE_ATOL, rtol=AFFINE_RTOL)
    # Non-symmetric matrix explicitly catches G transpose errors.
    assert affine["gradient"][0, 1] != affine["gradient"][1, 0]

from dataclasses import fields
import numpy as np
import pytest
from particle_3d.field import FlowSample


def assert_equivalent(field, positions):
    batch = field.sample_many(positions)
    for i, p in enumerate(positions):
        scalar = field.sample(p)
        for name in [f.name for f in fields(FlowSample)]:
            np.testing.assert_array_equal(getattr(batch, name)[i], getattr(scalar, name), strict=True)


def test_real_batch_matches_scalar_bitwise(real_field, random_real, faces):
    positions = np.vstack((random_real[1], real_field.points_m[[0, 100, 1000]],
                           [f["center"] for f in faces], [[0., 0., 0.], [1., 1., 1.]]))
    assert_equivalent(real_field, positions)


def test_synthetic_batch_matches_scalar_and_empty(affine):
    assert_equivalent(affine["field"], np.vstack((affine["positions"], [[20., 20., 20.]])))
    batch = affine["field"].sample_many(np.empty((0, 3)))
    assert batch.velocity_m_s.shape == (0, 3) and batch.velocity_gradient_s_inv.shape == (0, 3, 3)
    assert batch.tetra_id.dtype == np.int64 and batch.inside_lumen.dtype == bool


@pytest.mark.parametrize("positions", [[], [1, 2, 3], [[1., np.nan, 2]], np.zeros((2, 4))])
def test_invalid_batch_is_rejected(affine, positions):
    with pytest.raises(ValueError, match="positions_m"):
        affine["field"].sample_many(positions)

import numpy as np
import pytest


def test_expected_membership_and_all_invalid_values(classifications):
    for row in classifications:
        assert row["correct"], row
        if not row["expected_inside"]:
            assert row["tetra_id"] == -1 and row["all_physical_nan"], row


@pytest.mark.parametrize("position", [[np.nan, 0, 0], [np.inf, 0, 0], [1, 2], [[1, 2, 3]]])
def test_invalid_coordinate_is_explicit_error(affine, position):
    with pytest.raises(ValueError, match="position_m"):
        affine["field"].sample(position)


def test_outside_is_not_stationary_blood(affine):
    sample = affine["field"].sample([20., 20., 20.])
    assert not sample.inside_lumen and sample.tetra_id == -1
    assert np.isnan(sample.velocity_m_s).all() and np.isnan(sample.pressure_pa)
    assert np.isnan(sample.velocity_gradient_s_inv).all()
    assert np.isnan(sample.vorticity_s_inv).all() and np.isnan(sample.strain_rate_s_inv).all()

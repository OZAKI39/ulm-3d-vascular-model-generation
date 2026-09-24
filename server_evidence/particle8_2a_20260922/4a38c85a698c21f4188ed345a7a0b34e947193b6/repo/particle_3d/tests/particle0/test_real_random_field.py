import numpy as np


def test_all_real_fields_finite_and_match_constructed_point_weights(random_real, real_field):
    rows, positions, ids, weights, batch = random_real
    assert len(rows) == 400 and batch.inside_lumen.all()
    np.testing.assert_array_equal(batch.tetra_id, ids)
    expected_u = np.einsum("ni,nij->nj", weights, real_field.velocity_nodes_m_s[real_field.tetra[ids]])
    expected_p = np.einsum("ni,ni->n", weights, real_field.pressure_nodes_pa[real_field.tetra[ids]])
    # Per-cell error budget; reports retain actual errors separately.
    tau = real_field.geometry.weight_tolerance[ids]
    ubound = 8 * tau * np.max(np.abs(real_field.velocity_nodes_m_s[real_field.tetra[ids]]), axis=(1, 2))
    pbound = 8 * tau * np.max(np.abs(real_field.pressure_nodes_pa[real_field.tetra[ids]]), axis=1)
    assert np.all(np.max(np.abs(batch.velocity_m_s - expected_u), axis=1) <= ubound)
    assert np.all(np.abs(batch.pressure_pa - expected_p) <= pbound)
    for values in [batch.velocity_m_s, batch.pressure_pa, batch.velocity_gradient_s_inv, batch.vorticity_s_inv, batch.strain_rate_s_inv]:
        assert np.isfinite(values).all()
    np.testing.assert_array_equal(batch.strain_rate_s_inv, .5 * (batch.velocity_gradient_s_inv + batch.velocity_gradient_s_inv.transpose(0, 2, 1)))

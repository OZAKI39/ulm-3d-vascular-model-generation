"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.hydraulic_resistance import linear_radius_resistance
from network_1d0d.roi_mapping import insert_cuts, match_edge_position, clip_interval


def test_roi_exact_cut_insertion():
    xyz = np.array([[0., 0., 0.], [100e-6, 0., 0.]])
    radius = np.array([1e-6, 2e-6])
    result = insert_cuts([10, 20], xyz, radius, np.array([[0, 1]]), {0: [(.25, -2), (.7, -3)]})
    e = result['edges']; x = result['xyz_m']; r = result['radius_m']
    split_R = linear_radius_resistance(np.linalg.norm(x[e[:, 1]]-x[e[:, 0]], axis=1), r[e[:, 0]], r[e[:, 1]])
    assert sum(split_R) == pytest.approx(linear_radius_resistance(100e-6, 1e-6, 2e-6), rel=2e-14)
    np.testing.assert_allclose(x[2], [25e-6, 0., 0.], atol=1e-20)
    match = match_edge_position([0., 0., 0.], [100., 0., 0.], 1., 2., [25., 0., 0.], 1.25)
    assert match['fraction_parent_to_child'] == .25
    with pytest.raises(ValueError):
        match_edge_position([0., 0., 0.], [100., 0., 0.], 1., 2., [25., 1., 0.], 1.25)
    with pytest.raises(ValueError):
        insert_cuts([10, 20], xyz, radius, np.array([[0, 1]]), {0: [(.25, -2), (.25, -3)]})


def test_box_clip_both_endpoints_outside():
    assert clip_interval(np.array([-1., 0., 0.]), np.array([2., 0., 0.]), [0., -1., -1.], [1., 1., 1.]) == (1/3, 2/3)

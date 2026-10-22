"""Evidence sampling must not fabricate contact geometry or complete runs."""

import h5py
import numpy as np
from scripts.stage5a.compare import sample


def test_physical_clock_uses_actual_left_state_and_records_age(tmp_path):
    r = 1.78e-6
    x = np.array([[0.0, 0.0, 0.0], [2.01 * r, 0.0, 0.0], [0.0, 2.01 * r, 0.0]])
    p = tmp_path / "trajectory.h5"
    with h5py.File(p, "x") as h:
        h["time_s"] = [3.0, 3.00002]
        h["positions_m"] = np.array([x, x + [0.0, 1e-6, 0.0]])
        h["active"] = np.ones((2, 3), bool)
        h["radii_m"] = np.full(3, r)
        h["particle_id"] = [7, 13, 23]
    d = sample(p, np.array([3.0, 3.00001, 3.00002]))
    np.testing.assert_array_equal(d["positions_m"][0], d["positions_m"][1])
    np.testing.assert_array_equal(d["positions_m"][2], x + [0.0, 1e-6, 0.0])
    np.testing.assert_allclose(d["sample_age_s"], [0, 1e-5, 0], atol=1e-15)
    assert set(d["particle_ids"]) == {7, 13, 23}


def test_zero_common_clock_cannot_pass_scientific_comparison(tmp_path):
    from scripts.stage5a.compare import compare

    r = 1.78e-6
    x = np.array([[0.0, 0.0, 0.0], [2.01 * r, 0.0, 0.0], [0.0, 2.01 * r, 0.0]])
    files = []
    for name, end in [("old", 3.001), ("partial", 3.000002)]:
        p = tmp_path / (name + ".h5")
        files.append(p)
        with h5py.File(p, "x") as h:
            h["time_s"] = [3.0, end]
            h["positions_m"] = np.array([x, x + [0.0, 1e-8, 0.0]])
            h["active"] = np.ones((2, 3), bool)
            h["radii_m"] = np.full(3, r)
            h["particle_id"] = [0, 1, 2]
    report = compare(*files, tmp_path / "comparison", "partial")
    assert not report["pass_scientific"]
    assert report["compared_duration_s"] == 0
    assert report["Rg_endpoint_relative"] is None
    assert report["comparison_status"] == "NOT_EVALUABLE_NO_NONZERO_COMMON_CLOCK"

"""Unmatched legacy frames must not manufacture a continuous-path error bound."""

import h5py
import numpy as np
import pytest
from scripts.stage5a_r3.lineage import sample_all


def fixture(path, straight=False):
    with h5py.File(path, "w") as h:
        h["time_s"] = [3.0, 3.001]
        h["positions_m"] = [[[0, 0, 0]], [[2e-6, 0, 0]]]
        h["particle_id"] = [152]
        h["radii_m"] = [1.78e-6]
        h["active"] = [[True], [True]]
        if straight:
            h.attrs["actual_path"] = "piecewise linear accepted half steps; no coordinate interpolation"
            h.attrs["all_accepted_substeps_recorded"] = True


def test_legacy_chord_is_descriptive_at_unmatched_clock(tmp_path):
    p = tmp_path / "legacy.h5"
    fixture(p)
    data = sample_all(p, np.array([3.0, 3.0005, 3.001]))
    assert data["uncertainty_is_bound"].tolist() == [True, False, True]
    assert data["time"].tolist() == [3.0, 3.0, 3.001]
    assert data["positions"][1, 0, 0] == 0  # No interpolated midpoint.


def test_recorded_straight_segment_has_a_valid_hold_bound(tmp_path):
    p = tmp_path / "r3.h5"
    fixture(p, straight=True)
    data = sample_all(p, np.array([3.0005]))
    assert data["uncertainty_is_bound"].all()
    assert data["uncertainty_m"][0] >= 1e-6


def test_no_sampling_beyond_recorded_progress(tmp_path):
    p = tmp_path / "partial.h5"
    fixture(p)
    with pytest.raises(ValueError, match="actual stored interval"):
        sample_all(p, np.array([3.005]))

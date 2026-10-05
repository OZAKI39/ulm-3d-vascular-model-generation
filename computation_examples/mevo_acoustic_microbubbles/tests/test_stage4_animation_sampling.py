"""Saved snapshots and absorbing branch events must retain their own timing."""

import numpy as np
import pytest

# Rendering is an optional installation extra, including on compute-only hosts.
pytest.importorskip("imageio")
pytest.importorskip("PIL")
from scripts.stage4.render_motion import sample_case  # noqa: E402


@pytest.fixture
def crossing_case():
    # Nonconsecutive global IDs catch accidental use of an ID as an array index.
    return {
        "t": np.array([0., 1., 2., 3.]),
        "x": np.array([
            [[0., 0., 0.], [0., 0., 0.]],
            [[1., 2., 3.], [1., 0., 0.]],
            [[2., 4., 6.], [2., 0., 0.]],
            [[2., 4., 6.], [3., 0., 0.]],
        ]),
        "id_index": {71: 0, 900: 1},
        "events": [{"particle_id": 71, "crossing_time_s": 1.5,
                    "position_m": [2., 4., 6.]}],
    }


def test_event_between_snapshots_stops_at_recorded_crossing(crossing_case):
    sampled = sample_case(crossing_case, [1., 1.25, 1.5, 1.75, 3.])[0]
    np.testing.assert_array_equal(sampled[:, 0], [
        [1., 2., 3.], [1.5, 3., 4.5], [2., 4., 6.], [2., 4., 6.], [2., 4., 6.],
    ])
    # A surviving particle keeps its independent trajectory.
    np.testing.assert_array_equal(sampled[:, 1, 0], [1., 1.25, 1.5, 1.75, 3.])


def test_saved_states_and_event_on_snapshot_are_preserved(crossing_case):
    np.testing.assert_array_equal(sample_case(crossing_case, crossing_case["t"])[0], crossing_case["x"])
    crossing_case["events"][0]["crossing_time_s"] = 2.
    sampled = sample_case(crossing_case, [1.5, 2., 2.5])[0]
    np.testing.assert_array_equal(sampled[:, 0], [[1.5, 3., 4.5], [2., 4., 6.], [2., 4., 6.]])


@pytest.mark.parametrize("requested", [[-.001], [3.001], [float("nan")], [float("inf")]])
def test_display_cannot_extrapolate_or_use_invalid_time(crossing_case, requested):
    with pytest.raises(ValueError):
        sample_case(crossing_case, requested)

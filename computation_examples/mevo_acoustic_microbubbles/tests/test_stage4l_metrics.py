"""Gate statistics must distinguish rigid rotation from internal rearrangement."""

import json

import h5py
import numpy as np

from scripts.stage4l.analyze_dense import analyze, lifetime


def test_rigid_rotation_retains_graph_and_zero_rigid_residual(tmp_path):
    radius = 1.78e-6
    x = np.array([[0.0, 0.0, 0.0], [2.01 * radius, 0.0, 0.0], [1.005 * radius, 1.741 * radius, 0.0]])
    frames = []
    for angle in [0.0, 0.15, 0.3]:
        q = np.array(
            [[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]]
        )
        frames.append(x @ q + np.array([angle * 1e-5, 0.0, 0.0]))
    with h5py.File(tmp_path / "trajectories.h5", "x") as h:
        h.create_dataset("time_s", data=[0.0, 1e-5, 2e-5])
        h.create_dataset("positions_m", data=frames)
        h.create_dataset("active", data=np.ones((3, 3), bool))
        h.create_dataset("radii_m", data=np.full(3, radius))
        h.create_dataset("particle_id", data=[0, 1, 2])
    (tmp_path / "SUMMARY.json").write_text(json.dumps({"status": "COMPLETE"}))
    result = analyze(tmp_path, 5e-6)
    assert result["initial_cluster_size"] == 3
    assert result["sampling_mode"] == "last_accepted_state"
    assert result["max_sample_age_s"] > 0
    with np.load(result["metric_artifact"]) as sampled:
        np.testing.assert_allclose(
            sampled["radius_of_gyration_m"], sampled["radius_of_gyration_m"][0], rtol=1e-14
        )
    assert result["neighbor_turnover_mean"] == 0
    assert result["same_ID_centroid_removed_RMS_end_m"] > 1e-7
    assert result["rigid_fit_residual_end_m"] < 1e-20
    assert abs(result["gyration_relative_change"]) < 1e-14


def test_edge_lifetimes_censoring():
    a, b = (0, 1), (1, 2)
    result = lifetime([{a}, {a, b}, {b}], np.array([0.0, 0.1, 0.2]))
    assert result == [
        dict(pair=[0, 1], duration_s=0.2, left_censored=True, right_censored=False),
        dict(pair=[1, 2], duration_s=0.1, left_censored=False, right_censored=True),
    ]


def test_historical_primary_field_metadata():
    from scripts.stage4l.branch_runs import historical_field_path

    assert historical_field_path({"pzt_control": {"field": "B2_field.h5"}}, {}) == "B2_field.h5"
    assert (
        historical_field_path(
            {"effective_acoustic_field": {"path": "override.h5"}, "pzt_control": {"field": "B2_field.h5"}}, {}
        )
        == "override.h5"
    )


def test_native_graph_accounting_is_invariant_to_repeated_static_states():
    from scripts.stage4l.native_graph import graph_statistics

    before = np.array([[0.0, 0.0, 0.0], [0.9, 0.0, 0.0], [2.0, 0.0, 0.0]])
    after = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [2.0, 0.0, 0.0]])
    coarse = graph_statistics([before, after, after], [0.0, 0.1, 0.2], 1.0)
    fine = graph_statistics([before, before, after, after, after], [0.0, 0.05, 0.1, 0.15, 0.2], 1.0)
    assert coarse["edge_changes"] == fine["edge_changes"] == 2
    assert coarse["edge_changes_per_s"] == fine["edge_changes_per_s"] == 10
    assert coarse["observed_lifetime_mean_s"] == fine["observed_lifetime_mean_s"] == 0.1
    assert coarse["left_censored_lifetimes"] == fine["left_censored_lifetimes"] == 1
    assert coarse["right_censored_lifetimes"] == fine["right_censored_lifetimes"] == 1

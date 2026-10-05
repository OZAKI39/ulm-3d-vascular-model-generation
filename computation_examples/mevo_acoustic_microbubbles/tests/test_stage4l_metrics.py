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
    result = analyze(tmp_path, 1e-5)
    assert result["initial_cluster_size"] == 3
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

"""A failed or mismatched task must never produce a completion speedup."""

import h5py
import numpy as np
from scripts.stage5a_r3.gate import ratio, comparison_uses_run


def run(folder, wall, status="COMPLETE", dx=0):
    folder.mkdir()
    with h5py.File(folder / "trajectories.h5", "w") as h:
        h["time_s"] = [3.0, 3.005]
        h["positions_m"] = np.array([[[dx, 0, 0], [np.nan, np.nan, np.nan]]] * 2)
        h["active"] = [[True, False]] * 2
        h["particle_id"] = [152, 153]
        h["radii_m"] = [1.78e-6, 1.78e-6]
    return dict(
        run=str(folder),
        charged_wall_s=wall,
        summary=dict(status=status, completed_duration_s=0.005 if status == "COMPLETE" else 0.00431),
    )


def test_complete_matching_initial_state_uses_full_process_cost(tmp_path):
    old = run(tmp_path / "old", 500)
    new = run(tmp_path / "new", 100)
    r = ratio(old, new, 0.005)
    assert r["status"] == "PASS"
    assert r["speedup"] == 5
    assert r["same_initial"]


def test_no_extrapolation_from_failed_reference(tmp_path):
    old = run(tmp_path / "old", 500, "FAIL")
    new = run(tmp_path / "new", 100)
    r = ratio(old, new, 0.005)
    assert r["status"] == "NOT_EVALUABLE"
    assert r["speedup"] is None


def test_moved_initial_state_cannot_pass_timing_identity(tmp_path):
    old = run(tmp_path / "old", 500)
    new = run(tmp_path / "new", 100, dx=1e-15)
    r = ratio(old, new, 0.005)
    assert r["status"] == "NOT_EVALUABLE"
    assert r["speedup"] is None


def test_precision_evidence_is_bound_to_actual_run_not_same_case(tmp_path):
    comparison = dict(
        case="B2",
        left=str(tmp_path / "repeat0/trajectories.h5"),
        right=str(tmp_path / "refined/trajectories.h5"),
    )
    assert comparison_uses_run(comparison, dict(run=str(tmp_path / "repeat0")))
    assert comparison_uses_run(comparison, dict(run=str(tmp_path / "refined")))
    assert not comparison_uses_run(comparison, dict(run=str(tmp_path / "repeat1")))
    assert not comparison_uses_run(comparison, dict(run=None))

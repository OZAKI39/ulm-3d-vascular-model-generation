"""Small independent invariance and archive consistency checks, without GPU reruns."""

import numpy as np
import h5py
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.morphodynamics import rigid_fit, morphology_series, topology
from .analyze import BASE, BRANCHES, read, js
from .common import write


def main():
    x = np.array([[0, 0, 0], [4, 0, 0], [1, 3, 0], [1, 1, 2.0]]) * 1e-6
    angle = np.pi / 3
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    y = x @ rotation.T + np.array([10, 20, 30]) * 1e-6
    r, error = rigid_fit(x, y)
    assert error < 1e-19
    np.testing.assert_allclose(np.degrees(np.arccos((np.trace(r) - 1) / 2)), 60, atol=1e-10)
    changed = y.copy()
    changed[1, 0] += 1e-6
    assert rigid_fit(x, changed)[1] > 0.1e-6
    s = morphology_series(np.stack([x, y]), np.arange(4), np.array([[0, 0, 0, 0], [0, 0, 0, 1]]))
    assert s["dynamic_member_mask"][1].sum() == 3
    assert s["fixed_centroid_m"].shape == (2, 3)
    tetra = np.array([[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1.0]]) * 1e-6
    iso = morphology_series(tetra[None], np.arange(4), np.zeros((1, 4), dtype=int))
    assert not iso["principal_axis_identifiable"][0] and np.isnan(iso["principal_axis_angle_deg"][0])
    xx = np.array([[0, 0, 0], [3.6e-6, 0, 0], [10e-6, 0, 0]])
    t = np.array([0, 0.001, 0.0015, 0.004, 0.01])
    graphs = []
    for on in [True, False, True, False, False]:
        g = np.eye(3, dtype=bool)
        g[0, 1] = g[1, 0] = on
        graphs.append(connected_components(g, directed=False)[1])
    top, _, _ = topology(
        xx,
        t,
        np.array(graphs),
        np.array([0.001, 0.0015, 0.004]),
        np.array([[0, 0, 1, 0], [0, 0, 1, 1], [0, 0, 1, 0]]),
        0.01,
        0.002,
    )
    assert (
        top["counts"]["split_transient"] == 1
        and top["counts"]["split_persistent"] == 1
        and top["counts"]["merge_persistent"] == 1
    )
    state = np.load(BASE / "mature_v1/MATURE_STATE.npz")["positions_m"]
    for path in BRANCHES.values():
        d = read(path)
        np.testing.assert_array_equal(d["positions_m"][0], state)
        np.testing.assert_array_equal(d["particle_id"], np.arange(32))
        assert np.isfinite(d["positions_m"]).all()
        assert len(d["time_s"]) == 601 and abs(d["time_s"][-1] - 0.12) < 1e-14
        summary = js(path / "SUMMARY.json")
        assert summary["numerical_quality"] == "PASS"
        assert summary["macro_step_count"] == summary["actual_substep_count"] == 12000
        assert summary["quality"]["max_accepted_penetration_m"] <= 0.01 * 1.78e-6
    with h5py.File(BASE / "visualization/continuous_actual_states.h5") as h:
        d = read(BASE / "mature_v1/formation")
        e = read(BRANCHES["amplitude"])
        np.testing.assert_array_equal(h["positions_m"][:512], d["positions_m"][:512])
        np.testing.assert_array_equal(h["positions_m"][512:], e["positions_m"][1:])
        np.testing.assert_allclose(np.diff(h["time_s"][:]), 0.0002, rtol=0, atol=1e-14)
    write(
        BASE / "ANALYSIS_CHECKS.json",
        dict(
            status="PASS",
            checks=[
                "Rigid rotation/translation separated from deformation",
                "Dynamic vs fixed membership",
                "Degenerate principal axis masked",
                "Transient and persistent topology, including later invalidation",
                "All six branches start at identical saved state with exact IDs/time/quality",
                "Video source equals actual precursor and exact-state branch with uniform time",
            ],
            scope="Offline semantic and data consistency checks; not physical validation",
        ),
    )
    print("PASS: six offline invariance/data checks")


if __name__ == "__main__":
    main()

"""Rigid fitting and group membership changes must preserve every sphere."""

import numpy as np
from scripts.stage4l import magnified_geometry as g


class SphericalWall:
    def __init__(self, config=None):
        pass

    def bind_mesh_surface(self, mesh):
        pass

    def wall_gap_and_normal(self, x, r):
        distance = np.linalg.norm(x, axis=1)
        normal = x / np.maximum(distance[:, None], 1e-30)
        return 0.001 - distance - r, normal, np.ones(len(x), bool)


def test_wall_fit_moves_whole_group_without_deforming_it():
    positions = np.array([[0.001, 0, 0], [0.001, 0.0002, 0], [0.001, -0.0002, 0]])
    radii = np.array([4e-5, 5e-5, 3e-5])
    original = positions.copy()
    shift, _ = g.fit_rigid_group(SphericalWall(), positions, radii, np.zeros(3))
    result = positions + shift
    assert SphericalWall().wall_gap_and_normal(result, radii)[0].min() >= 0
    np.testing.assert_allclose(result - result[0], original - original[0], atol=1e-18, rtol=0)
    np.testing.assert_array_equal(positions, original)
    assert np.linalg.norm(shift) > 0


def test_merge_retains_all_ids_and_same_scale_for_radii_and_vectors(monkeypatch):
    monkeypatch.setattr(g, "BifurcationGeometry", SphericalWall)
    monkeypatch.setattr(g.pv, "read", lambda path: None)
    x = np.array(
        [
            [[0.0008, 0, 0], [0.000805, 0, 0], [0.00081, 0.00001, 0]],
            [[0.000801, 0, 0], [0.000804, 0, 0], [0.000807, 0.000003, 0]],
        ]
    )
    before = x.copy()
    d = dict(
        config={"geometry": {}},
        surface="unused",
        x=x,
        t=np.array([0, 0.02]),
        ids=np.array([7, 19, 41]),
        radii=np.array([1e-6, 1.5e-6, 1e-6]),
        id_index={7: 0, 19: 1, 41: 2},
        states=[dict(active_mask=np.ones(3, bool)) for _ in range(2)],
        clusters=[[dict(particle_ids=[7, 19]), dict(particle_ids=[41])], [dict(particle_ids=[7, 19, 41])]],
    )
    samples, groups, report = g.prepare(d, 24)
    np.testing.assert_array_equal(x, before)
    np.testing.assert_array_equal(samples["radii_m"], d["radii"] * 24)
    separated = np.linalg.norm(samples["positions_m"][0, :2] - samples["positions_m"][0, 2], axis=1)
    gaps = separated - samples["radii_m"][:2] - samples["radii_m"][2]
    assert gaps.min() >= 5e-6 - 1e-12
    for k, rows in enumerate(groups):
        assert sorted(pid for row in rows for pid in row["particle_ids"]) == [7, 19, 41]
        for row in rows:
            ix = [d["id_index"][pid] for pid in row["particle_ids"]]
            native = x[k, ix]
            display = samples["positions_m"][k, ix]
            np.testing.assert_allclose(display - display[0], 24 * (native - native[0]), atol=1e-18, rtol=0)
    assert report["min_display_surface_gap_m"] >= 0

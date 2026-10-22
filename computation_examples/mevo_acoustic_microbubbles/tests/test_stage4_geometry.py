import numpy as np
import pytest
import yaml
from mevo_mb.coupling.bifurcation import BifurcationGeometry, BranchClassifier, Branch


def geometry():
    with open("configs/stage4/bifurcation.yaml") as stream:
        return BifurcationGeometry(yaml.safe_load(stream)["geometry"])


def test_geometry_connected_symmetric_and_finite():
    g = geometry()
    parent = np.column_stack([np.linspace(-.013, 0, 80), np.zeros(80), np.zeros(80)])
    for direction in g.directions.values():
        path = np.vstack([parent, np.linspace(0, .011, 80)[:, None] * direction])
        assert g.contains(path, 1.78e-6).all()
    rng = np.random.default_rng(28)
    x = rng.uniform([-.013, -.006, -.0014], [.011, .006, .0014], (700, 3))
    mirrored = x.copy()
    mirrored[:, 1] *= -1
    assert np.max(np.abs(g.wall_distance(x) - g.wall_distance(mirrored))) < 1e-14
    assert np.isfinite(g.wall_distance(x)).all()
    assert not g.contains(np.array([[-.02, 0, 0], [.02, 0, 0], [0, 0, .01]])).any()


def test_classifier_first_crossing_and_permanent_particle_identity():
    g = geometry()
    c = BranchClassifier(g, [73, 19, 88])
    a, b = g.directions.values()
    old = np.array([a*.009, b*.009, a*.011])
    new = np.array([a*.011, b*.011, a*.009])
    velocity = (new-old)/.2
    events = c.observe([73, 19, 88], old, new, velocity, velocity, 1., 1.2)
    assert len(events) == 2
    assert c.assignment == {73: Branch.BRANCH_A, 19: Branch.BRANCH_B, 88: Branch.UNASSIGNED}
    assert abs(events[0]["crossing_time_s"] - 1.1000001) < 1e-12
    assert np.allclose(events[0]["position_m"], a*(.010+1e-9), atol=1e-15)
    assert not c.observe([73], (b*.009)[None], (b*.011)[None], b[None], b[None], 2, 3)
    assert c.assignment[73] == Branch.BRANCH_A


def test_classifier_rejects_missed_disk_and_endpoint_guess():
    g = geometry()
    c = BranchClassifier(g, [0])
    a = g.directions[Branch.BRANCH_A]
    old, new = a*.009+[0,0,.002], a*.011+[0,0,.002]
    assert not c.observe([0], old[None], new[None], a[None], a[None], 0, 1)
    assert c.summary()["N_unassigned"] == 1
    with pytest.raises(ValueError):
        BranchClassifier(g, [4, 4])


def test_projection_and_geometric_wall_slip():
    g = geometry()
    x = np.array([[-.006, .0015, 0.]])
    radius = np.array([1.78e-6])
    feasible, hit = g.project(x, radius)
    assert hit.all() and g.contains(feasible, radius-1e-12).all()
    velocity = g.constrain_velocity(feasible, radius, np.array([[.003, .004, 0.]]))
    assert velocity[0, 0] == pytest.approx(.003)
    assert abs(velocity[0, 1]) < 1e-9


def test_partial_experiment_preserves_accepted_state(tmp_path):
    """A query failure must preserve evidence instead of publishing a complete trajectory."""
    import h5py
    import json
    from mevo_mb.bubbles.population import BubblePopulation
    from mevo_mb.coupling.bifurcation_experiment import BifurcationExperimentTracker
    from mevo_mb.io.artifacts import new_run
    from mevo_mb.utils.config import load_config

    config = load_config("configs/demo.yaml")
    config["geometry"] = dict(geometry().config,kind="bifurcation")
    config["particles"].update(primary_enabled=False,secondary_enabled=False,duration_s=.001,dt_s=.0001)
    config["experiment_scope"] = "SOFTWARE_TEST_ONLY, force-free, deliberate interpolation failure"

    class Field:
        frequency_Hz = 2.5e6

        def interpolate_complex_pressure(self,x):
            return np.zeros(len(x),complex)

        def interpolate_complex_gradient(self,x):
            return np.zeros_like(x,dtype=complex)

    class BrokenFlow:
        def __init__(self):
            self.calls = 0

        def velocity(self,x,t):
            self.calls += 1
            if self.calls > 1:
                raise ValueError("Deliberate rejected interpolation query")
            return np.tile([.001,0.,0.],(len(x),1))

    population = BubblePopulation(np.array([[-.002,.0002,0.],[-.002,-.0002,0.]]),
                                  np.full(2,1.78e-6),np.arange(2))
    out = new_run(tmp_path,"partial",config)
    with pytest.raises(ValueError,match="Deliberate"):
        BifurcationExperimentTracker(config,Field(),BrokenFlow()).run(population,out)
    with h5py.File(out / "particles/trajectories.h5") as h:
        assert not h.attrs["complete"]
        np.testing.assert_array_equal(h["positions_m"][0],population.positions_m)
    assert json.loads((out / "EXECUTION_FAILURE.json").read_text())["status"] == "INCOMPLETE"
    with np.load(out / "particles/last_accepted_state.npz") as saved:
        np.testing.assert_array_equal(saved["positions_m"],population.positions_m)


def test_mesh_distance_cache_handles_mutated_queries_and_results():
    import pyvista as pv
    g = geometry()
    g.bind_mesh_surface(pv.Sphere(radius=.001).triangulate())
    x = np.array([[.0002,0.,0.]])
    distance,gradient = g.level_and_gradient(x)
    expected = distance.copy()
    distance[:] = 999.
    gradient[:] = 999.
    repeated,_ = g.level_and_gradient(x)
    np.testing.assert_array_equal(repeated,expected)
    assert g.distance_cache_hits == 1
    x[0,0] = .0003
    moved,_ = g.level_and_gradient(x)
    assert moved[0] > expected[0]
    assert g.distance_cache_misses == 2


@pytest.mark.parametrize("per_particle", [False, True])
def test_wall_certificate_skips_only_proven_noop_and_preserves_exact_queries(per_particle):
    import pyvista as pv
    from scripts.stage4.geometry_certificate import install

    original, certified = geometry(), geometry()
    surface = pv.Sphere(radius=.001).triangulate()
    for g in [original, certified]:
        g.bind_mesh_surface(surface)
    counts = install(certified, per_particle=per_particle)
    radii = np.full(2, 1.78e-6)
    warm = np.array([[.0002, 0., 0.], [0., .0003, 0.]])
    for g in [original, certified]:
        g.wall_distance(warm)
    x = warm + np.array([[1e-6, 2e-6, 0.], [0., 0., 3e-6]])
    a, b = original.project(x, radii), certified.project(x, radii)
    for actual, expected in zip(b, a, strict=True):
        np.testing.assert_array_equal(actual, expected)
    assert counts["certified_inactive"] == 1
    # A bound must never be passed off as an exact distance/gradient cache entry.
    np.testing.assert_array_equal(certified._distance_cache[0], warm)
    np.testing.assert_array_equal(certified.wall_distance(x), original.wall_distance(x))
    # A possible wall collision must still use the original correction.
    outside = x.copy()
    outside[0] = [.0011, 0., 0.]
    a, b = original.project(outside, radii), certified.project(outside, radii)
    for actual, expected in zip(b, a, strict=True):
        np.testing.assert_array_equal(actual, expected)
    assert b[1][0] and counts["certified_inactive"] == 1
    assert counts["partial_batches"] == int(per_particle)


@pytest.mark.parametrize("per_particle", [False, True])
def test_wall_certificate_falls_back_without_compatible_euclidean_cache(per_particle):
    import pyvista as pv
    from scripts.stage4.geometry_certificate import install

    g = geometry()
    g.bind_mesh_surface(pv.Sphere(radius=.001).triangulate())
    counts = install(g, per_particle=per_particle)
    x, r = np.array([[.0002, 0., 0.]]), np.array([1.78e-6])
    g.project(x, r)  # No prior cache.
    assert counts["certified_inactive"] == 0
    g._distance_cache[2][:] *= .5  # Cannot certify unit-gradient Euclidean distance.
    g.project(x+1e-6, r)
    assert counts["certified_inactive"] == 0
    g.project(np.empty((0, 3)), np.empty(0))
    assert counts["certified_inactive"] == 0

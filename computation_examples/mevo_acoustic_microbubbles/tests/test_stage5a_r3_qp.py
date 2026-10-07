import numpy as np
import pytest
from scipy import sparse
from mevo_mb.bubbles.motion_qp_r2 import certificate
from mevo_mb.bubbles.motion_qp_r3 import NativeQPR3
from mevo_mb.bubbles.motion_r3 import MotionR3, SemiImplicitR3
from mevo_mb.bubbles.motion_backend_r3 import MotionR3Optimized


def test_wrong_or_stale_free_hint_is_not_a_certificate():
    a = sparse.diags([1.0, 2.0, 3.0], format="csc")
    g = sparse.eye(3, format="csc")
    q = NativeQPR3()
    for b in [np.array([-1.0, 3.0, -2.0]), np.array([4.0, -1.0, 6.0])]:
        u, lam, _ = q.solve(a, b, g, np.zeros(3), free_hint=np.zeros(3))
        assert certificate(a, b, g, np.zeros(3), u, lam)["pass_gate"]
        assert np.allclose(u, np.maximum(b / a.diagonal(), 0), atol=1e-9)


@pytest.mark.parametrize("workers,components", [(1, True), (3, True), (1, False)])
def test_changing_state_component_global_and_parallel_agreement(workers, components):
    n = 12
    r = np.full(n, 1.78e-6)
    x = np.zeros((n, 3))
    x[:, 0] = np.arange(n) * 10e-6
    x[1::2, 0] = x[::2, 0] + 2.01 * r[0]
    ids = np.arange(n)
    wall = (np.full(n, 1e-3), np.tile([0.0, 1.0, 0.0], (n, 1)), np.ones(n, bool))
    old = MotionR3(0.001)
    new = MotionR3Optimized(0.001, components=components, component_workers=workers)
    accepted = new.cache
    for k in range(4):
        force = np.zeros_like(x)
        force[::2, 0] = 1e-9 * (1 + 0.1 * k)
        force[1::2, 0] = -force[::2, 0]
        a = old.solve(x, r, force, np.zeros_like(x), wall)
        from mevo_mb.bubbles.motion_r2 import anchored_constraints

        c = anchored_constraints(x, x, r, wall, wall, 1e-6, ids, 1e-11, a.velocity)
        u0, _, d0 = old.constrained_stage(a, c, ids, "changing")
        u, cache, d = new.constrained_stage(a, c, ids, "changing")
        assert new.cache is accepted
        assert d["kkt"] <= 1e-9 and d0["kkt"] <= 1e-9
        assert np.max(abs(u - u0)) < 1e-9
        x = x + 1e-7 * u
    if new.executor:
        new.executor.shutdown()


def test_finite_wall_feature_rows_do_not_add_lubrication_terms():
    # Two finite perpendicular faces at equal contact distance. Resistance is
    # still the single supplied nearest-WALL normal, whereas G covers both faces.

    class Sweep:
        triangles = np.array([[[0, 0, 0], [0, 4, 0], [0, 0, 4]], [[0, 0, 0], [4, 0, 0], [0, 0, 4]]], float)

        def candidates(self, p, q, reach):
            return np.array([0, 1])

    b = MotionR3(0.001)
    x = np.array([[1.0, 1.0, 1.0]])
    r = np.array([1.0])
    wall = (np.zeros(1), np.array([[-1.0, 0.0, 0.0]]), np.ones(1, bool))
    step = SemiImplicitR3(b, None, None, Sweep())
    before = b.solve(x, r, np.zeros_like(x), np.ones_like(x), wall)
    c = step.feature_constraints(x, r, wall, 0.001, np.array([152]), before.velocity)
    after = b.solve(x, r, np.zeros_like(x), np.ones_like(x), wall)
    assert np.array_equal(before.matrix.toarray(), after.matrix.toarray())
    assert sum(k[0] == "wall_feature" for k in c.ids) == 2
    assert np.allclose(c.matrix.toarray()[-2:], [[1, 0, 0], [0, 1, 0]])

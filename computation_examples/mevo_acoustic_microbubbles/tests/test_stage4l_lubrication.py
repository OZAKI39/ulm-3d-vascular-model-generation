"""Analytic and limiting-case checks for the scoped normal-resistance model."""

import numpy as np
import pytest
from scipy.sparse.linalg import cg

from mevo_mb.bubbles.translation_backends import (
    LegacyStokesTranslation,
    LubricationCoupledTranslation,
)
from mevo_mb.coupling.geometry import Cylinder

R, MU = 1.78e-6, 0.001


def pair(eps):
    return np.array([[0.0, 0.0, 0.0], [(2 + eps) * R, 0.0, 0.0]])


def wall(n=2, gap=np.inf):
    return np.full(n, gap), np.tile([0.0, 0.0, 1.0], (n, 1)), np.ones(n, bool)


def solve(x, force=None, flow=None, w=None):
    n = len(x)
    model = LubricationCoupledTranslation(MU)
    result = model.solve(
        x,
        np.full(n, R),
        np.zeros((n, 3)) if force is None else force,
        np.zeros((n, 3)) if flow is None else flow,
        wall(n) if w is None else w,
    )
    return model, result


def test_T1_single_farfield_bitwise():
    x, r = np.zeros((1, 3)), np.array([R])
    f, u = np.array([[1e-11, -2e-12, 4e-11]]), np.array([[0.003, 0.002, -0.001]])
    legacy = LegacyStokesTranslation(MU).solve(x, r, f, u)
    _, new = solve(x, f, u)
    np.testing.assert_array_equal(new.velocity, legacy.velocity)
    assert new.diagnostics["solver"] == "exact_legacy_bypass"


def test_T2_outside_activation_zero():
    model = LubricationCoupledTranslation(MU)
    np.testing.assert_array_equal(model.coefficient(np.array([0.1, 0.11, 1]) * R, R), 0)
    for eps in [0.1, 0.11, 1.0]:
        x, r, u = pair(eps), np.full(2, R), np.arange(6.0).reshape(2, 3)
        new = model.solve(x, r, u * 1e-12, u)
        old = LegacyStokesTranslation(MU).solve(x, r, u * 1e-12, u)
        np.testing.assert_allclose(new.velocity, old.velocity, rtol=5e-15, atol=1e-16)


@pytest.mark.parametrize("eps", [0.1, 0.05, 0.01, 0.001])
def test_T3_pair_analytic_squeeze_factor_two(eps):
    model, result = solve(pair(eps), force=np.array([[1e-11, 0, 0], [-1e-11, 0, 0]]))

    def phi(e):
        return 1 / (2 * e) - 9 / 20 * np.log(e) - 3 / 56 * e * np.log(e)

    expected = 3 * np.pi * MU * R * (phi(eps) - phi(0.1))
    assert model.coefficient(eps * R, R) == pytest.approx(expected, abs=1e-22)
    drag = 6 * np.pi * MU * R
    expected_v = 1e-11 / (drag + 2 * expected)
    np.testing.assert_allclose(result.velocity[:, 0], [expected_v, -expected_v], rtol=1e-12)
    if eps == 0.001:
        # Leading pair coefficient uses relative speed U_i-U_j, not individual speed.
        leading = 3 * np.pi * MU * R**2 / (2 * eps * R)
        assert expected / leading == pytest.approx(1.0, rel=0.02)


def test_T4_symmetry_reciprocity_spd_and_cg():
    x = np.array([[0, 0, 0], [2.01 * R, 0, 0], [R, 1.75 * R, 0]])
    model, result = solve(x, force=np.arange(9.0).reshape(3, 3) * 1e-12)
    a = result.matrix.toarray()
    np.testing.assert_allclose(a, a.T, atol=1e-22)
    assert np.linalg.eigvalsh(a).min() > 0
    b = np.arange(9.0) * 1e-12
    vc, status = cg(result.matrix / result.scale, b / result.scale, rtol=1e-13, atol=0)
    assert status == 0
    np.testing.assert_allclose(vc, result.velocity.ravel(), rtol=1e-11, atol=1e-14)
    correction = a - np.diag(np.repeat(result.drag, 3))
    force = (correction @ np.arange(9.0)).reshape(3, 3)
    np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-20)


def test_T5_common_translation_no_pair_force():
    v = np.tile([0.01, -0.02, 0.003], (2, 1))
    _, result = solve(pair(0.001), flow=v)
    np.testing.assert_allclose(result.velocity, v, rtol=1e-12)
    np.testing.assert_allclose(
        (result.matrix @ v.ravel()).reshape(2, 3), result.drag[:, None] * v, atol=1e-22
    )


@pytest.mark.parametrize("eps", [0.1, 0.05, 0.01, 0.001])
def test_T6_T7_wall_normal_and_tangent(eps):
    force = np.array([[1e-11, 2e-11, 3e-11]])
    model, result = solve(np.zeros((1, 3)), force=force, w=wall(1, eps * R))

    def phi(e):
        return 1 / e - 1 / 5 * np.log(e) - 1 / 21 * e * np.log(e)

    aw = 6 * np.pi * MU * R * (phi(eps) - phi(0.1))
    drag = 6 * np.pi * MU * R
    np.testing.assert_allclose(result.velocity, force / [drag, drag, drag + aw], rtol=1e-12)
    assert model.coefficient(eps * R, R, "wall") == pytest.approx(aw, abs=1e-22)
    if eps == 0.001:
        assert aw / (6 * np.pi * MU * R**2 / (eps * R)) == pytest.approx(1, rel=0.02)


@pytest.mark.parametrize("kind", ["pair", "wall"])
def test_T8_activation_continuity(kind):
    m = LubricationCoupledTranslation(MU)
    assert m.coefficient(0.1 * R, R, kind) == 0
    assert m.coefficient((0.1 - 1e-10) * R, R, kind) / (6 * np.pi * MU * R) < 2e-8
    assert m.coefficient((0.1 + 1e-10) * R, R, kind) == 0


@pytest.mark.parametrize("kind", ["pair", "wall"])
def test_T9_floor_saturation_not_equilibrium(kind):
    m = LubricationCoupledTranslation(MU)
    a = m.coefficient(np.array([0.001, 0.0005, 0, -0.001]) * R, R, kind)
    assert np.isfinite(a).all()
    np.testing.assert_array_equal(a, np.full(4, a[0]))
    _, result = solve(pair(0.0005), force=np.array([[1e-11, 0, 0], [-1e-11, 0, 0]]))
    assert result.velocity[0, 0] > 0
    safe, _ = m.contact(result, pair(0.0005), np.full(2, R), result.velocity, 1e-11, horizon=0.001)
    assert safe[0, 0] > 0  # positive gap may close, never glued at floor
    assert (safe[0, 0] - safe[1, 0]) * 0.001 <= 0.0005 * R + 1e-18


@pytest.mark.parametrize("mode", ["closing", "separating", "tangential", "subfloor"])
def test_T10_pair_timestep(mode):
    m = LubricationCoupledTranslation(MU)
    eps = 0.0005 if mode == "subfloor" else 0.05
    v = np.array([[1.0, 0, 0], [-1.0, 0, 0]])
    if mode == "separating":
        v = -v
    if mode == "tangential":
        v = v[:, [2, 0, 1]]
    dt = m.timestep(pair(eps), np.full(2, R), v, wall(), 0.001)
    assert dt > 0
    assert dt == pytest.approx(0.5 * eps * R / 2 if mode == "closing" else 0.001)


@pytest.mark.parametrize("mode", ["closing", "separating", "tangential", "subfloor"])
def test_T10_wall_timestep(mode):
    m = LubricationCoupledTranslation(MU)
    v = np.array([[0.0, 0.0, 1.0]])
    if mode == "separating":
        v = -v
    if mode == "tangential":
        v = v[:, [2, 0, 1]]
    eps = 0.0005 if mode == "subfloor" else 0.05
    dt = m.timestep(np.zeros((1, 3)), np.array([R]), v, wall(1, eps * R), 0.001)
    assert dt > 0
    assert dt == pytest.approx(0.5 * eps * R if mode == "closing" else 0.001)


def test_predictive_pair_outside_activation():
    m = LubricationCoupledTranslation(MU)
    x = pair(10.0)
    v = np.array([[1.0, 0, 0], [-1.0, 0, 0]])
    dt = m.timestep(x, np.full(2, R), v, wall(), 0.001)
    assert dt == pytest.approx(0.5 * 10 * R / 2)
    assert np.linalg.norm((x + dt * v)[1] - (x + dt * v)[0]) - 2 * R > 0.1 * R
    # Common translation does not alter the gap limiter.
    assert m.timestep(x, np.full(2, R), v + 1000, wall(), 0.001) == pytest.approx(dt)


def test_unequal_radii_rejected():
    with pytest.raises(NotImplementedError, match="UNEQUAL_RADIUS"):
        LubricationCoupledTranslation(MU).solve(
            pair(0.05), np.array([R, R * 1.01]), np.zeros((2, 3)), np.zeros((2, 3))
        )


def test_A_metric_contact_and_wall():
    x = pair(0)
    force = np.array([[1e-11, 2e-12, 2e-11], [-1e-11, 0, 0]])
    model, result = solve(x, force=force, w=wall(2, 0))
    safe, diag = model.contact(result, x, np.full(2, R), result.velocity, 1e-11)
    assert safe[1, 0] - safe[0, 0] >= -1e-12
    assert np.max(safe[:, 2]) <= 1e-12
    assert diag["contact_kkt_residual"] < 1e-9
    assert safe[0, 1] == pytest.approx(force[0, 1] / (6 * np.pi * MU * R))
    # Orthogonality / objective against feasible perturbations.
    grad = result.matrix @ (safe - result.velocity).ravel()
    for v in [np.zeros((2, 3)), safe + np.array([[-0.1, 0, -0.1], [0.1, 0, -0.1]])]:
        assert grad @ (v - safe).ravel() >= -1e-22


def test_geometry_normal_gap():
    geometry = Cylinder(0.001, 0.01)
    x = np.array([[0, 0.001 - R - 0.01 * R, 0], [0.005 - R, 0, 0]])
    g, n, valid = geometry.wall_gap_and_normal(x, np.full(2, R))
    np.testing.assert_allclose(g, [0.01 * R, 0], atol=1e-18)
    np.testing.assert_array_equal(n, [[0, 1, 0], [1, 0, 0]])
    assert valid.all()


def test_lubrication_independent_of_acoustic_switches():
    from mevo_mb.bubbles.translation_backends import make_translation_backend

    backend = make_translation_backend(
        dict(
            medium={"viscosity_Pa_s": MU},
            particles={
                "translation_backend": "lubrication_coupled",
                "primary_enabled": False,
                "secondary_enabled": False,
            },
        )
    )
    flow = np.array([[0.01, 0.0, 0.0], [-0.01, 0.0, 0.0]])
    result = backend.solve(pair(0.01), np.full(2, R), np.zeros((2, 3)), flow, wall())
    assert result.diagnostics["pair_active"] == 1
    assert 0 < result.velocity[0, 0] < flow[0, 0]


def test_vtk_wall_gap_reuses_exact_cache():
    import pyvista as pv
    from mevo_mb.coupling.bifurcation import BifurcationGeometry

    geometry = BifurcationGeometry(
        dict(
            parent_length_m=0.014,
            parent_radius_m=0.0013,
            daughter_length_m=0.012,
            daughter_radius_m=0.001,
            blending_length_m=0.0002,
            branch_angle_deg=30,
            classification_distance_m=0.01,
            bifurcation_plane_x_m=0.0,
            classification_tolerance_m=1e-9,
        )
    )
    geometry.bind_mesh_surface(pv.Sphere(radius=0.001, theta_resolution=20, phi_resolution=20))
    x = np.array([[0.0, 0.0, 0.0009], [0.0, 0.0008, 0.0]])
    value, gradient = geometry.level_and_gradient(x)
    misses = geometry.distance_cache_misses
    gap, normal, valid = geometry.wall_gap_and_normal(x, np.full(2, R))
    assert geometry.distance_cache_misses == misses
    np.testing.assert_array_equal(gap, -value / np.linalg.norm(gradient, axis=1) - R)
    np.testing.assert_array_equal(normal, gradient / np.linalg.norm(gradient, axis=1)[:, None])
    assert valid.all()


@pytest.mark.parametrize(
    "fixture",
    [
        "stage4l_contact_qp",
        "stage4l_contact_qp_b2",
        "stage4l_contact_qp_secondary",
        "stage4l_contact_qp_cycle",
    ],
)
def test_observed_nearly_dependent_contact_qp_qr_refinement(fixture):
    """Replay the rejected B3 primary RK stage without weakening KKT tolerance."""
    from pathlib import Path
    from scipy.sparse import load_npz
    from scipy.sparse.linalg import splu
    from mevo_mb.bubbles.translation_backends import TranslationResult

    folder = Path(__file__).parent / "data" / fixture
    with np.load(folder / "contact_input.npz") as data:
        d = {key: data[key] for key in data.files}
    a = load_npz(folder / "matrix.npz")
    scale = float(d["scale"])
    result = TranslationResult(
        d["velocity_m_s"],
        d["velocity_m_s"],
        a,
        d["drag"],
        None,
        None,
        None,
        (d["wall_gap_m"], d["wall_normal"], d["wall_valid"]),
        factor=splu(a / scale),
        scale=scale,
    )
    model = LubricationCoupledTranslation(MU)
    velocity, diag = model.contact(
        result,
        d["positions_m"],
        d["radii_m"],
        d["velocity_m_s"],
        float(d["tolerance_m"]),
        None if np.isnan(float(d["horizon_s"])) else float(d["horizon_s"]),
    )
    assert diag["contact_kkt_residual"] < 1e-10
    assert diag["contact_primal_stationarity"] < 1e-10
    assert diag["contact_refinement"] in (
        "whitened_primal_QR_active_set",
        "feasible_primal_QR_active_set",
    )
    delta = (velocity - d["velocity_m_s"]).ravel()
    zero = -d["velocity_m_s"].ravel()
    assert delta @ (a @ delta) <= zero @ (a @ zero) * (1 + 1e-10)


def test_touching_pair_tangential_motion_is_not_locked():
    free = np.array([[0.0, 0.01, 0.0], [0.0, -0.01, 0.0]])
    model, result = solve(pair(0.0), flow=free)
    velocity, diag = model.contact(result, pair(0.0), np.full(2, R), result.velocity, 1e-11, horizon=0.001)
    np.testing.assert_allclose(velocity, free, rtol=1e-13, atol=1e-15)


def test_returned_resistance_components_obey_dimensional_force_balance():
    force = np.array([[1e-11, 2e-12, 3e-11], [-2e-11, 4e-12, 1e-11]])
    model, result = solve(pair(0.01), force=force, w=wall(2, 0.01))
    pair_force, wall_force = result.lubrication_forces()
    np.testing.assert_allclose(pair_force.sum(axis=0), 0, atol=1e-25)
    np.testing.assert_allclose(wall_force[:, :2], 0, atol=1e-25)
    np.testing.assert_allclose(
        force + pair_force + wall_force,
        result.drag[:, None] * result.velocity,
        rtol=1e-11,
        atol=1e-24,
    )
    assert result.diagnostics["pair_subfloor"] == 0
    assert result.diagnostics["wall_subfloor"] == 0
    tangent = np.array([[0, 0.2, 0], [0, -0.2, 0]])
    pp, pw = result.lubrication_forces(tangent)
    np.testing.assert_array_equal(pp, np.zeros((2, 3)))
    np.testing.assert_array_equal(pw, np.zeros((2, 3)))


def test_nonfinite_force_is_rejected_with_json_safe_failure_evidence():
    import json

    model = LubricationCoupledTranslation(MU)
    force = np.full((2, 3), np.nan)
    with pytest.raises(RuntimeError, match="LINEAR_SOLVER_RESIDUAL_FAILURE"):
        model.solve(pair(0.01), np.full(2, R), force, np.zeros((2, 3)), wall())
    assert model.last_failure["relative_residual"] is None
    assert "relative_residual" in model.last_failure["nonfinite_diagnostics"]
    json.dumps(model.last_failure, allow_nan=False)

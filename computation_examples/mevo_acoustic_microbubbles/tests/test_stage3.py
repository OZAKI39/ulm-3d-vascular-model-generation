import ast
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid

from mevo_mb.bubbles.contact import contact_candidates, project_contacts, constrain_pair_velocities
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from mevo_mb.bubbles.validity import PairPhysicsRegime
from mevo_mb.coupling.geometry import Cylinder

spec = importlib.util.spec_from_file_location(
    "stage3_history", Path("validation/history_force_sensitivity/model.py")
)
history = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = history
spec.loader.exec_module(history)


def test_secondary_direction_scaling_reciprocity_zero():
    model = ClassicalFarFieldBjerknes()
    r = np.ones(2) * 1e-6
    distances = np.array([10, 20, 30, 50, 100]) * 1e-6
    forces = []
    for d in distances:
        x = np.array([[0.0, 0, 0], [d, 0, 0]])
        same = model.forces(x, r, np.array([1 + 1j, 1 + 1j]) * 1e-19, 2e6, 1000)
        opposite = model.forces(x, r, np.array([1 + 1j, -1 - 1j]) * 1e-19, 2e6, 1000)
        np.testing.assert_allclose(same[0], -same[1], rtol=1e-14)
        assert same[0, 0] > 0 and opposite[0, 0] < 0
        assert not model.forces(x, r, np.array([1e-19, 0]), 2e6, 1000).any()
        forces.append(same[0, 0])
    assert np.polyfit(np.log(distances), np.log(forces), 1)[0] == pytest.approx(-2, abs=1e-12)
    scale = np.array([1, 2, 3, 4])
    f = [model.forces(x, r, np.ones(2) * 1e-19 * p, 2e6, 1000)[0, 0] for p in scale]
    assert np.polyfit(np.log(scale), np.log(f), 1)[0] == pytest.approx(2, abs=1e-12)


def test_warning_policy_does_not_regularize_force():
    r = np.full(2, 1e-6)
    v = np.ones(2) * 1e-19
    m = ClassicalFarFieldBjerknes(enforce_far_field=False)
    f = []
    for d in [2e-6, 3e-6, 10e-6]:
        x = np.array([[0.0, 0, 0], [d, 0, 0]])
        f.append(m.forces(x, r, v, 1e6, 1000)[0, 0] * d * d)
    np.testing.assert_allclose(f, f[0], rtol=1e-14)
    with pytest.raises(ValueError):
        m.forces(np.zeros((2, 3)), r, v, 1e6, 1000)
    regime = PairPhysicsRegime()
    for d, key in [(2e-6, "contact"), (3e-6, "near_field"), (6e-6, "transition"), (12e-6, "far_field")]:
        assert regime.population(np.array([[0.0, 0, 0], [d, 0, 0]]), r)[key + "_fraction"] == 1


@pytest.mark.parametrize("seed", range(5))
def test_kdtree_exact_candidates_and_resolution(seed):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-15e-6, 15e-6, (60, 3))
    r = rng.uniform(0.5e-6, 1.5e-6, 60)
    for gap in [-1e-11, 0, 1e-11, 2e-6]:
        np.testing.assert_array_equal(
            contact_candidates(x, r, gap, "brute"), contact_candidates(x, r, gap, "ckdtree")
        )
    cylinder = Cylinder(radius_m=0.0015, length_m=0.02)
    a = project_contacts(x, r, cylinder, 1e-11, 300, "brute")
    b = project_contacts(x, r, cylinder, 1e-11, 300, "ckdtree")
    np.testing.assert_array_equal(a[0], b[0])
    assert a[2:] == b[2:]
    velocity = rng.normal(size=x.shape) * 0.01
    np.testing.assert_array_equal(
        constrain_pair_velocities(a[0], r, velocity, 1e-11, 300, "brute"),
        constrain_pair_velocities(b[0], r, velocity, 1e-11, 300, "ckdtree"),
    )


def test_variable_radius_history_integrates_linear_q_exactly():
    t = np.linspace(0, 0.001, 301)
    r = 1e-6 * (1 + 0.1 * np.sin(t / 0.001 * 2 * np.pi))
    theta = cumulative_trapezoid(1 / r**2, t, initial=0)
    q = 3e-12 * theta
    np.testing.assert_allclose(
        history.history_integral(theta, q), 6e-12 * np.sqrt(theta), rtol=1e-13, atol=1e-25
    )


def test_history_fixed_radius_added_mass_and_manufactured_solution():
    r = 1e-6
    rho = 1000
    mu = 0.001
    t = np.linspace(0, 5e-6, 2001)
    force = np.full(len(t), 1e-11)
    force[0] = 0
    tau = (0.5 * rho * 4 * np.pi * r**3 / 3) / (6 * np.pi * mu * r)
    b = history.solve(t, np.full(len(t), r), force)
    exact = 1e-11 / (6 * np.pi * mu * r) * (1 - np.exp(-t / tau))
    assert np.max(abs(b["velocity_m_s"] - exact)) / max(exact) < 0.001
    # Analytic Abel integral of q=a*theta: 2*a*sqrt(theta).
    theta = t / r**2
    a = 1e-14
    q = a * theta
    v = q / r
    mass = 0.5 * rho * 4 * np.pi * r**3 / 3
    forcing = (
        mass * a / r**3
        + 6 * np.pi * mu * r * v
        + 6 * np.pi * mu / np.sqrt(np.pi * mu / rho) * 2 * a * np.sqrt(theta)
    )
    c = history.solve(t, np.full(len(t), r), forcing, history=True)
    np.testing.assert_allclose(c["velocity_m_s"], v, rtol=1e-10, atol=1e-13)
    assert not history.solve(t, np.full(len(t), r), np.zeros(len(t)), history=True)["velocity_m_s"].any()


def test_production_has_no_history_dependency_and_no_old_gates():
    for path in Path("src/mevo_mb").rglob("*.py"):
        if "/validation/" in str(path):
            continue
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert "history_force_sensitivity" not in (node.module or "")
            if isinstance(node, ast.Import):
                assert all("history_force_sensitivity" not in v.name for v in node.names)
    from mevo_mb.validation.gates import REQUIRES

    assert not any(
        "ZHANG" in name or "PAIR_REFERENCE" in name for names in REQUIRES.values() for name in names
    )


def test_completed_stage3_evidence():
    import json

    root = Path("results/stage3")
    if not (root / "cylinder/STATE.json").exists():
        pytest.skip("Formal cylinder acceptance evidence not generated yet")
    for kind in ["secondary", "clustering", "flow", "history", "cylinder"]:
        data = json.loads((root / kind / "STATE.json").read_text())
        assert data["status"] == "PASS", kind
    data = json.loads((root / "cylinder/STATE.json").read_text())
    assert data["tests"]["phase_changes_distribution"] and data["tests"]["flow_only_analytic"]


def test_contact_tolerance_boundary_preserves_original_strict_criterion():
    r = np.full(2, 1e-6)
    tolerance = 1e-11
    threshold = r.sum() - tolerance
    cylinder = Cylinder(radius_m=0.0015, length_m=0.02)
    for d in [np.nextafter(threshold, 0), threshold, np.nextafter(threshold, np.inf)]:
        x = np.array([[0.0, 0, 0], [d, 0, 0]])
        amount = r.sum() - d
        expected = x.copy()
        if amount > tolerance:
            expected[0, 0] -= amount / 2
            expected[1, 0] += amount / 2
        for backend in ["brute", "ckdtree"]:
            actual = project_contacts(x, r, cylinder, tolerance, 10, backend)[0]
            np.testing.assert_array_equal(actual, expected)
    with pytest.raises(ValueError):
        ClassicalFarFieldBjerknes(min_distance_over_sum=np.nan)

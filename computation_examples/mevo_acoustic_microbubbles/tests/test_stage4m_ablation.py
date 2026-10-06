"""M1--M10: fixed-scale sensitivity policy, including the frozen Stage 4L path."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from mevo_mb.bubbles.secondary_ablation import near_contact_secondary_ablation as weight

R = 1.78e-6


def backend(policy, **kw):
    return ClassicalFarFieldBjerknes(enforce_far_field=False, near_contact_secondary_ablation=policy, **kw)


def pair_force(policy, gap, phase=0):
    x = np.array([[0.0, 0.0, 0.0], [2 * R + gap, 0.0, 0.0]])
    return backend(policy).forces(x, np.full(2, R), np.array([1e-18, 1e-18 * np.exp(1j * phase)]), 1e6, 1000)


def test_M1_full_above_activation():
    np.testing.assert_array_equal(weight(np.array([0.1, 0.11, 1, 100]) * R, R, "A2"), 1)


def test_M2_zero_below_floor():
    np.testing.assert_array_equal(weight(np.array([-0.1, 0, 0.0009, 0.001]) * R, R, "A2"), 0)


def test_M3_continuity():
    for endpoint, value in [(0.001, 0), (0.1, 1)]:
        np.testing.assert_allclose(
            weight((endpoint + np.array([-1e-7, 0, 1e-7])) * R, R, "A2"), value, atol=2e-14
        )


@pytest.mark.parametrize("endpoint", [0.001, 0.1])
def test_M4_first_derivative(endpoint):
    # Finite differences of the implemented function, derivative in h/R.
    delta = 1e-5
    values = weight((endpoint + np.array([-delta, delta])) * R, R, "A2")
    assert abs((values[1] - values[0]) / (2 * delta)) < 6e-7


@pytest.mark.parametrize("endpoint", [0.001, 0.1])
def test_M5_second_derivative(endpoint):
    delta = 1e-6
    values = weight((endpoint + np.array([-delta, 0, delta])) * R, R, "A2")
    assert abs((values[2] - 2 * values[1] + values[0]) / delta**2) < 0.02


@pytest.mark.parametrize("policy", ["A0", "A1", "A2"])
@pytest.mark.parametrize("phase", [0, 2.5])
def test_M6_pair_reciprocity(policy, phase):
    f = pair_force(policy, 0.03 * R, phase)
    np.testing.assert_array_equal(f[0], -f[1])


@pytest.mark.parametrize("harmonics", [1, 5])
@pytest.mark.parametrize("diagnostics", [False, True])
def test_M7_frozen_A0_bitwise(harmonics, diagnostics):
    path = Path(__file__).parent / "data/stage4m/stage4l_forces_frozen.py"
    spec = importlib.util.spec_from_file_location("stage4l_frozen", path)
    frozen = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = frozen
    spec.loader.exec_module(frozen)
    rng = np.random.default_rng(20261006)
    x = rng.normal(size=(31, 3)) * 5e-6
    v = (rng.normal(size=(31, harmonics)) + 1j * rng.normal(size=(31, harmonics))) * 1e-18
    if harmonics == 1:
        v = v[:, 0]
    args = (x, np.full(31, R), v, 1e6, 1000)
    expected = frozen.ClassicalFarFieldBjerknes(enforce_far_field=False).forces(*args)
    np.testing.assert_array_equal(backend("A0", record_pair_diagnostics=diagnostics).forces(*args), expected)


def test_M8_subfloor_only():
    for gap in [-0.0001, 0, 0.0009]:
        np.testing.assert_array_equal(pair_force("A1", gap * R), 0)
    for gap in [0.0011, 0.02, 0.1, 3]:
        np.testing.assert_array_equal(pair_force("A1", gap * R), pair_force("A0", gap * R))
    np.testing.assert_array_equal(weight(np.array([0.001]) * R, R, "A1"), 0)


def test_M9_A2_only_inside_activation():
    for gap in [0.1001, 1, 3]:
        np.testing.assert_array_equal(pair_force("A2", gap * R), pair_force("A0", gap * R))
    np.testing.assert_allclose(
        pair_force("A2", 0.0505 * R), 0.5 * pair_force("A0", 0.0505 * R), rtol=1e-13, atol=0
    )


def test_M10_secondary_disabled_identical():
    # Real tracker evaluate path, manufactured radial/field fixture, no scientific rerun.
    import yaml
    from mevo_mb.bubbles.tracker import ParticleTracker

    class Field:
        frequency_Hz = 1e6

        def interpolate_complex_pressure(self, x):
            return np.full(len(x), 1000 + 0j)

        def interpolate_complex_gradient(self, x):
            return np.tile([1e5 + 0j, 2e5 + 0j, 0j], (len(x), 1))

    class Flow:
        def velocity(self, x, t):
            return np.tile([0.001, 0.0, 0.0], (len(x), 1))

    states = []
    for policy in ["A0", "A1", "A2"]:
        config = yaml.safe_load(Path("configs/demo.yaml").read_text())
        config["particles"].update(
            secondary_enabled=False,
            translation_backend="lubrication_coupled",
            near_contact_secondary_ablation=policy,
        )
        tracker = ParticleTracker(config, Field(), Flow())
        tracker.pair = backend(policy)
        tracker.pair.forces = lambda *args: (_ for _ in ()).throw(
            AssertionError("Disabled secondary evaluated")
        )
        x = np.array([[0.0, 0.0, 0.0], [2.03 * R, 0.0, 0.0]])
        states.append(tracker.evaluate(x, np.full(2, R), 0))
    for state in states[1:]:
        for key in state:
            np.testing.assert_array_equal(state[key], states[0][key])


def test_diagnostics_pair_norm_not_net_and_null():
    m = backend("A2", record_pair_diagnostics=True)
    x = np.array([[0.0, 0.0, 0.0], [2.0505 * R, 0.0, 0.0]])
    f = m.forces(x, np.full(2, R), np.full(2, 1e-18 + 0j), 1e6, 1000)
    d = m.pair_diagnostics()
    assert d["effective_classical_ratio"] == pytest.approx(0.5)
    assert d["S_effective_N"] == pytest.approx(np.linalg.norm(f[0]))
    assert d["unordered_pair_count"] == d["taper_zone_pair_count"] == 1
    m.forces(x, np.full(2, R), np.zeros(2, complex), 1e6, 1000)
    assert m.pair_diagnostics()["effective_classical_ratio"] is None


def test_pairwise_weighting_not_particle_net_scalar():
    x = np.array([[0.0, 0.0, 0.0], [2.03 * R, 0.0, 0.0], [0.0, 3 * R, 0.0]])
    v, radii = np.full(3, 1e-18 + 0j), np.full(3, R)
    expected = np.zeros_like(x)
    for i in range(3):
        for j in range(i + 1, 3):
            f = backend("A2").forces(x[[i, j]], radii[[i, j]], v[[i, j]], 1e6, 1000)
            expected[[i, j]] += f
    np.testing.assert_allclose(backend("A2").forces(x, radii, v, 1e6, 1000), expected, rtol=1e-15, atol=1e-28)


@pytest.mark.gpu
@pytest.mark.parametrize("policy", ["A0", "A1", "A2"])
def test_numpy_cupy_pair_policy_consistency(policy):
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() < 1:
        pytest.skip("No CUDA device")
    x = np.array([[0.0, 0.0, 0.0], [2.0005 * R, 0.0, 0.0], [0.0, 2.04 * R, 0.0], [0.0, 0.0, 3.0 * R]])
    v = np.array([[1.0, 0.2], [0.8j, 0.1], [-1.0, 0.3j], [0.4 + 0.8j, 0.2j]]) * 1e-18
    args = (x, np.full(4, R), v, 1e6, 1000)
    gpu = ClassicalFarFieldBjerknes(
        "cupy", enforce_far_field=False, near_contact_secondary_ablation=policy, record_pair_diagnostics=True
    )
    cpu = backend(policy, record_pair_diagnostics=True)
    np.testing.assert_allclose(gpu.forces(*args).get(), cpu.forces(*args), rtol=2e-14, atol=1e-27)
    g, c = gpu.pair_diagnostics(), cpu.pair_diagnostics()
    for key in ["full_weight_pair_count", "taper_zone_pair_count", "zeroed_pair_count"]:
        assert g[key] == c[key]
    assert g["S_effective_N"] == pytest.approx(c["S_effective_N"], rel=2e-14, abs=1e-27)


@pytest.mark.parametrize(
    "ratios,stable,expected",
    [
        ((5, 1, 1, 5), True, True),
        ((4.99, 4.99, 5, 1), True, False),
        ((5, 1, 1, 4.99), True, False),
        ((10, 10, 10, 10), False, False),
    ],
)
def test_engineering_gate_requires_each_case_and_stability(tmp_path, monkeypatch, ratios, stable, expected):
    from scripts.stage4m import analyze

    monkeypatch.setattr(analyze, "BASE", tmp_path)
    (tmp_path / "comparisons").mkdir()
    rows = {}
    for i, case in enumerate(["B2", "B3"]):
        rows[case] = {
            "A0": {"summary": dict(accepted_steps=10000, wall_s=1000)},
            "A2": {
                "summary": dict(
                    accepted_steps=10000 / ratios[2 * i],
                    wall_s=1000 / ratios[2 * i + 1],
                    status="COMPLETE" if stable else "INCOMPLETE_RUNTIME",
                    completed_duration_s=0.005,
                    error=None,
                    projection={"multi_radius_corrections": 0},
                    source_inputs_unchanged=True,
                ),
                "all_positions_finite": True,
                "all_steps_positive": True,
                "max_rg_ratio": 1.01,
                "translation": {"failure": None},
            },
        }
    assert analyze.dense_gate(rows)["continue_formation"] is expected


def test_empty_active_population_has_no_stale_pair_diagnostics():
    from types import SimpleNamespace
    from scripts.stage4m.common import recorded_pair_diagnostics

    tracker = SimpleNamespace(config={"medium": {"density_kg_m3": 1000}})
    d = recorded_pair_diagnostics(tracker, 0)
    assert d["unordered_pair_count"] == 0
    assert d["S_classical_N"] == d["S_effective_N"] == 0
    assert d["effective_classical_ratio"] is None
    assert d["full_weight_pair_fraction"] is None

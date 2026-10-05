from dataclasses import replace
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("independent_zhang", Path("references/zhang2011_minimal.py"))
ref = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ref
spec.loader.exec_module(ref)


def test_independent_residual_matches_random_state_not_shared_algebra():
    from mevo_mb.validation.zhang2011 import Zhang2011EncapsulatedPairReference, Zhang2011Parameters

    p = ref.Parameters(D0=50e-6)
    m = Zhang2011EncapsulatedPairReference(Zhang2011Parameters())
    rng = np.random.default_rng(25)
    for _ in range(10):
        y = np.array(
            [
                p.R10 * (1 + rng.uniform(-0.1, 0.1)),
                rng.normal(),
                0.0,
                rng.normal() * 0.01,
                p.R20 * (1 + rng.uniform(-0.1, 0.1)),
                rng.normal(),
                p.D0,
                rng.normal() * 0.01,
            ]
        )
        t = rng.random() / p.f
        acc = ref.acceleration(t, y, p)
        yp = y[[0, 4, 2, 6, 1, 5, 3, 7]]
        np.testing.assert_allclose(acc, m.rhs(t, yp)[4:], rtol=3e-12, atol=1e-7)
        assert np.max(abs(ref.residual(t, y, acc, p))) < 1e-10


def test_secondary_time_origin_and_coordinate_conventions():
    from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes

    n = 4096
    t = np.arange(n) / n
    h = np.arange(1, 4)
    v = np.array([[1 + 2j, 0.2 + 0.1j, -0.1j], [-0.5 + 2j, 0.3 - 0.2j, 0.2]]) * 1e-18

    def avg(phase):
        dv = np.real(np.exp(2j * np.pi * t[:, None] * h + 1j * phase * h) @ (v * (2j * np.pi * h)).T)
        return np.mean(dv[:, 0] * dv[:, 1])

    np.testing.assert_allclose(avg(0), avg(np.pi / 2), rtol=1e-13, atol=0)
    x = np.array([[0.0, 0, 0], [30e-6, 0, 0]])
    r = np.full(2, 1e-6)
    f = ClassicalFarFieldBjerknes().forces(x, r, v, 1e6, 1000)
    reflected = ClassicalFarFieldBjerknes().forces(-x, r, v, 1e6, 1000)
    np.testing.assert_allclose(reflected, -f, rtol=1e-13, atol=0)
    assert np.dot(f[0], x[1] - x[0]) > 0 and np.dot(reflected[0], -x[1] + x[0]) > 0


@pytest.mark.parametrize("kwargs", [{"Pa": np.nan}, {"R10": -1e-6}, {"D0": 1e-6}, {"chi": -1}, {"f": 0}])
def test_invalid_reference_input(kwargs):
    with pytest.raises(ValueError):
        replace(ref.Parameters(), **kwargs)


@pytest.mark.deprecated_reference
def test_completed_forensic_matrix_preserves_failure():
    import json

    path = Path("results/stage2_5/zhang_forensic/STATE.json")
    if not path.exists():
        pytest.skip("Run forensic evidence fixture first")
    data = json.loads(path.read_text())
    assert data["independent_implementation_agrees"]
    cases = {r["case"]: r for r in data["cases"]}
    assert cases["A_attraction"]["rows"][0]["calculated_sign"] == 1
    assert cases["B_repulsion"]["rows"][0]["calculated_sign"] == -1
    assert cases["C_fig3_failed"]["rows"][0]["result"] == "FAIL_SIGN"
    assert data["status"] == "UNRESOLVED" and not data["fig3_sign_pass"]


def test_synthetic_vtu_space_time_si_units(tmp_path):
    import pyvista as pv
    from mevo_mb.flow.fields import VTKFlowField

    points = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1.0, 0], [0.0, 0, 1.0]])
    files = []
    for i in [0, 1]:
        mesh = pv.UnstructuredGrid(np.array([4, 0, 1, 2, 3]), np.array([pv.CellType.TETRA]), points.copy())
        mesh["v"] = points * (i + 1)
        path = tmp_path / f"{i}.vtu"
        mesh.save(path)
        files.append(path)
    flow = VTKFlowField(
        files, [0.0, 2.0], velocity_array="v", coordinate_scale_to_m=0.001, velocity_scale_to_m_s=0.01
    )
    np.testing.assert_allclose(
        flow.velocity([[0.0001, 0.0002, 0.0003]], 1.0), [[0.0015, 0.003, 0.0045]], rtol=1e-12, atol=1e-15
    )
    with pytest.raises(ValueError):
        flow.velocity([[0.002, 0, 0]], 1.0)


def test_completed_three_model_far_limit_evidence():
    import json

    path = Path("results/stage2_5/pair_asymptotics/STATE.json")
    if not path.exists():
        pytest.skip("Run pair asymptotic benchmark first")
    data = json.loads(path.read_text())
    assert data["status"] == "COMPARISON_COMPLETE"
    for name in ["zhang", "doinikov"]:
        rows = sorted(
            [
                r
                for r in data["rows"]
                if r["model"] == name
                and r["pressure_Pa"] == 100
                and r["frequency_Hz"] == 300000
                and r["sound_speed_m_s"] == 1e9
            ],
            key=lambda r: r["ratio"],
        )
        assert len(rows) == 5 and all(r["calculated_sign"] == r["classical_sign"] for r in rows)
        errors = [r["coefficient_relative_error"] for r in rows]
        assert all(a > b for a, b in zip(errors, errors[1:])) and errors[-1] < 0.001
        assert rows[-1]["hydrodynamic_relative_error"] < 0.001
    # Formal c→infinity convergence does not certify finite-c / coated production.
    from mevo_mb.validation.regime_maps import PairPhysicsRegimeMap

    assert (
        PairPhysicsRegimeMap(path).classify(1.5e-6, 2e-6, 320e-6, 2.95e6, 100, "uncoated")
        == "OUTSIDE_VALIDATED_RANGE"
    )


def test_real_adaptive_midpoint_and_holdout_evidence():
    import json

    path = Path("results/stage2_5/response_table/STATE.json")
    if not path.exists():
        pytest.skip("Run adaptive database benchmark first")
    data = json.loads(path.read_text())
    run = Path(data["run"])
    tree = json.loads((run / "CHECKPOINT.json").read_text())
    assert data["status"] == "PASS_EMPIRICAL_DOMAIN" and data["coverage_fraction"] == pytest.approx(1)
    assert data["resume_check"]["new_apecss_calls"] == 0
    assert all(c["status"] == "VALID_CELL" and max(c["errors"].values()) <= 0.05 for c in tree["leaves"])
    hold = json.loads((run / "HOLDOUTS.json").read_text())
    assert hold["seed"] == 20261004
    assert len(hold["independent_uniform_queries"]) == 256
    assert all(
        r["status"] == "PASS" and max(r["errors"].values()) <= 0.05
        for r in hold["independent_uniform_queries"]
    )

from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from mevo_mb.bubbles.apecss import ApecssBubbleResponseModel, RadialCase, complex_harmonics
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution, load_parameter_cases
from mevo_mb.bubbles.response_table import ResponseTableBubbleModel, ResponseTableBoundsError

PARAMS = Path("configs/stage2/sonovue_parameters.json")


def test_parameter_identity_units_provenance():
    cases = load_parameter_cases(PARAMS)
    assert cases["br14_table1"].agent == "BR14"
    p = cases["sonovue_table2"]
    assert p.agent == "SonoVue" and p.radius_m.value == 1.78e-6
    assert p.shell_viscosity_kg_s.value == 3.2e-9
    assert p.shell_elasticity_N_m.source_doi == "10.1121/1.3242346"
    with pytest.raises(ValueError, match="N/m"):
        replace(p, shell_elasticity_N_m=replace(p.shell_elasticity_N_m, unit="MPa"))


def test_radial_damping_does_not_invalidate_acoustic_cache():
    from mevo_mb.acoustics.kwave_backend import acoustic_identity
    from mevo_mb.utils.config import load_config

    c = load_config("configs/demo.yaml")
    before = acoustic_identity(c)
    c["medium"]["radial_effective_viscosity_Pa_s"] = 0.002
    assert acoustic_identity(c) == before
    c["medium"]["sound_speed_m_s"] += 1
    assert acoustic_identity(c) != before


def test_radius_interpolation_bounds_and_reproducibility():
    dist = SonoVueParameterDistribution(PARAMS)
    r = np.array([1.0, 1.78, 2.5, 3.0]) * 1e-6
    a = dist.evaluate_batch(r)
    assert a["elasticity"].shape == (4,)
    assert np.all(a["viscosity"] > 0)
    np.testing.assert_array_equal(
        dist.evaluate_batch(r, uncertainty_draw=True, seed=29)["viscosity"],
        dist.evaluate_batch(r, uncertainty_draw=True, seed=29)["viscosity"],
    )
    with pytest.raises(ValueError, match="Outside"):
        dist.evaluate_batch([0.75e-6, 4e-6])
    assert "digitization" in dist.parameter_set(1.78e-6).shell_viscosity_kg_s.uncertainty_kind


def test_phasor_extraction_dc_harmonics_and_time_origin():
    f = 2.5e6
    t = (9.37 + np.arange(8 * 512) / 512) / f
    expected = np.array([2.3 + 1.7j, -0.2 + 0.7j, 0.3 - 0.1j])
    y = 13 + np.real(np.exp(2j * np.pi * f * t[:, None] * np.arange(1, 4)) @ expected)
    np.testing.assert_allclose(complex_harmonics(t, y, f, 3), expected, rtol=1e-12, atol=1e-13)
    with pytest.raises(ValueError, match="integer"):
        complex_harmonics(t[:-1], y[:-1], f)


@pytest.fixture
def affine_table(tmp_path):
    path = tmp_path / "response.h5"
    axes = (np.array([1e-6, 2e-6, 3e-6]), np.array([1e6, 2e6]), np.array([1e4, 2e4]))
    r, f, p = np.meshgrid(*axes, indexing="ij")
    v = 1e-18 * (r / 1e-6 + 2 * f / 1e6 + 3j * p / 1e4)
    with h5py.File(path, "w") as h5:
        h5.attrs.update(
            schema_version=1,
            complete=True,
            metadata_json=json.dumps(
                {"medium": {"density_kg_m3": 1000.0, "sound_speed_m_s": 1500.0, "viscosity_Pa_s": 0.002}}
            ),
        )
        for name, a in zip(("radius_m", "frequency_Hz", "pressure_Pa"), axes):
            h5[name] = a
        h5["harmonics_m3"] = np.stack((v, 0.1 * v), axis=-1)
        h5["node_valid"] = np.ones(r.shape, bool)
        h5["harmonic_resolved"] = np.ones(r.shape, bool)
        h5["cell_valid"] = np.ones((2, 1, 1), bool)
        h5["max_R_over_R0"] = np.full(r.shape, 1.1)
        h5["min_R_over_R0"] = np.full(r.shape, 0.9)
    return path


def test_response_table_affine_batch_phase_and_rejection(affine_table):
    model = ResponseTableBubbleModel(affine_table)
    r = np.array([1e-6, 1.8e-6, 3e-6])
    f, p = 1.6e6, np.array([1e4, 1.7e4, 2e4])
    expected = 1e-18 * (r / 1e-6 + 2 * f / 1e6 + 3j * p / 1e4)
    np.testing.assert_allclose(model.evaluate_batch(r, f, p), expected, rtol=1e-14, atol=0)
    np.testing.assert_allclose(
        model.response(r, f, p * np.exp(0.4j), {}), expected * np.exp(0.4j), rtol=1e-14, atol=0
    )
    harmonic = model.harmonic_response(r, f, p * np.exp(0.4j))["harmonics_m3"]
    np.testing.assert_allclose(harmonic[:, 1], 0.1 * expected * np.exp(0.8j), rtol=1e-14, atol=0)
    for bad in (0, 2.1e4, np.nan):
        with pytest.raises(ResponseTableBoundsError):
            model.evaluate_batch(r, f, bad)
    with pytest.raises(ValueError, match="medium mismatch"):
        model.response(r, f, p, {"density_kg_m3": 1060.0})
    np.testing.assert_allclose(
        model.response(r, f, p, {"viscosity_Pa_s": 0.001, "radial_effective_viscosity_Pa_s": 0.002}),
        expected,
        rtol=1e-14,
        atol=0,
    )


def test_table_invalid_cells_and_nodes_are_not_silently_used(affine_table):
    with h5py.File(affine_table, "r+") as h5:
        h5["cell_valid"][:] = False
    model = ResponseTableBubbleModel(affine_table)
    with pytest.raises(ResponseTableBoundsError, match="Unvalidated"):
        model.evaluate_batch(1.5e-6, 1.5e6, 1.5e4)
    assert np.isfinite(model.evaluate_batch(1e-6, 1e6, 1e4))
    with h5py.File(affine_table, "r+") as h5:
        h5["node_valid"][0, 0, 0] = False
    with pytest.raises(ResponseTableBoundsError, match="Nonperiodic"):
        ResponseTableBubbleModel(affine_table).evaluate_batch(1e-6, 1e6, 1e4)


@pytest.mark.apecss
def test_real_apecss_wrapper(tmp_path):
    exe = Path("third_party/build/apecss/ultrasound_apecss")
    if not exe.exists():
        pytest.skip("APECSS binary not built in this environment")
    model = ApecssBubbleResponseModel(exe)
    case = RadialCase(1.78e-6, 2.5e6, 5e3, 0.3, 3.2e-9, 0.036)
    meta, h, cycles = model.solve(case, tmp_path / "real")
    assert meta["returncode"] == 0 and meta["metrics"]["periodic"]
    assert h.shape == (8,) and np.all(np.isfinite(h))
    # The actual upstream excitation's pressure phasor is +i*pa.
    raw = np.loadtxt(next((tmp_path / "real").glob("RPAR_*.txt")))
    np.testing.assert_allclose(
        raw[:, 7],
        case.ambient_pressure_Pa - case.pressure_Pa * np.sin(2 * np.pi * case.frequency_Hz * raw[:, 1]),
        rtol=1e-12,
    )
    extracted = complex_harmonics(cycles["time_s"], cycles["volume_m3"], case.frequency_Hz)
    np.testing.assert_allclose(h, extracted / (1j ** np.arange(1, 9)), rtol=1e-12, atol=0)
    with pytest.raises(RuntimeError, match="offline"):
        model.response(1e-6, 1e6, 1e4, {})


def test_strict_tracker_table_and_preserved_invalid_state(affine_table, tmp_path):
    from mevo_mb.acoustics.basis import ComplexAcousticField
    from mevo_mb.bubbles.population import BubblePopulation
    from mevo_mb.bubbles.tracker import ParticleTracker
    from mevo_mb.bubbles.validity import PhysicsValidityError
    from mevo_mb.flow.fields import ZeroFlowField
    from mevo_mb.io.artifacts import new_run
    from mevo_mb.utils.config import load_config

    c = load_config("configs/demo.yaml")
    c["bubble_response"] = {"model": "response_table", "path": str(affine_table)}
    c["medium"].update(density_kg_m3=1000.0, sound_speed_m_s=1500.0, viscosity_Pa_s=0.002)
    c["particles"].update(near_field_policy="strict", duration_s=1e-7, dt_s=5e-8, secondary_enabled=True)
    axes = tuple(np.linspace(-0.01, 0.01, 3) for _ in range(3))
    field = ComplexAcousticField(axes, np.full((3, 3, 3), 1.5e4 + 0j), np.zeros((3, 3, 3, 3), complex), 1.5e6)
    for distance, valid in ((25e-6, True), (3e-6, False)):
        pop = BubblePopulation(
            np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]), np.full(2, 1e-6), np.arange(2)
        )
        run = new_run(tmp_path, "strict", c)
        tracker = ParticleTracker(c, field, ZeroFlowField())
        if valid:
            summary = tracker.run(pop, run)
            assert summary["physics_validity"]["physics_status"] == "partially validated"
            assert summary["physics_validity"]["within_response_table"]
        else:
            with pytest.raises(PhysicsValidityError):
                tracker.run(pop, run)
            assert (run / "particles/rejected_state.npz").exists()
            assert json.loads((run / "particles/SUMMARY.json").read_text())["status"] == "OUTSIDE_VALIDITY"
            with h5py.File(run / "particles/trajectories.h5") as h:
                assert not h.attrs.get("complete", False)

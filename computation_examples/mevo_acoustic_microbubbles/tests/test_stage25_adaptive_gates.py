from dataclasses import replace
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from mevo_mb.bubbles.adaptive_table import (
    ResponseNodeCache,
    node_key,
    AdaptiveResponseTableBubbleModel,
    CORNERS,
)
from mevo_mb.bubbles.apecss import RadialCase
from mevo_mb.validation.gates import ValidationGateManager, ValidationGateError, validation_watermark
from mevo_mb.validation.regime_maps import PairPhysicsRegimeMap


def case():
    return RadialCase(1.78e-6, 2.5e6, 5e3, 0.3, 3.2e-9, 0.036)


def test_cache_failed_node_no_repeat_and_crash_recovery(tmp_path, monkeypatch):
    from mevo_mb.bubbles.apecss import ApecssBubbleResponseModel

    exe = tmp_path / "apecss"
    exe.write_bytes(b"fixture")
    cache = ResponseNodeCache(tmp_path / "cache", exe)
    calls = []

    def fail(self, c, out):
        calls.append(1)
        raise RuntimeError("injected APECSS termination")

    monkeypatch.setattr(ApecssBubbleResponseModel, "solve", fail)
    a = cache.get(case())
    b = cache.get(case())
    assert a["status"] == b["status"] == "FAILED" and len(calls) == 1 and b["cache_hit"]
    c = replace(case(), pressure_Pa=6000)
    folder = cache.root / node_key(c, cache.binary_hash)
    attempt = folder / "attempt_0000"
    attempt.mkdir(parents=True)
    np.savez(attempt / "cycles.npz", harmonics_cosine_m3=np.ones(8) * 1e-18, radius_m=np.array([c.radius_m]))
    (attempt / "RESULT.json").write_text(
        json.dumps(dict(status="COMPLETE", metrics={"valid": True}, returncode=0, elapsed_s=1))
    )
    recovered = cache.get(c)
    assert recovered["recovered_completed_attempt"] and len(calls) == 1
    (attempt / "cycles.npz").unlink()
    with pytest.raises(ValueError, match="Incomplete"):
        cache.get(c)


def test_key_captures_numerics_medium_shell_and_provenance_is_not_physics():
    a = case()
    key = node_key(a, "binary")
    assert node_key(replace(a, pressure_Pa=5000), "binary") == key
    for name, val in [
        ("cycles", 96),
        ("samples_per_cycle", 1024),
        ("elasticity_N_m", 0.31),
        ("density_kg_m3", 1010),
        ("tolerance", 1e-11),
    ]:
        assert node_key(replace(a, **{name: val}), "binary") != key
    assert node_key(replace(a, provenance="same model, alternate citation formatting"), "binary") == key


@pytest.fixture
def adaptive_fixture(tmp_path):
    path = tmp_path / "adaptive.h5"
    lo = np.array([1e-6, 1e6, 1e4])
    hi = np.array([2e-6, 2e6, 2e4])
    pts = lo + CORNERS * (hi - lo)
    v = (pts[:, 0] / 1e-6 + 1j * pts[:, 1] / 1e6 + pts[:, 2] / 1e4) * 1e-18
    with h5py.File(path, "w") as h:
        h.attrs.update(
            schema_version=2,
            complete=True,
            metadata_json=json.dumps(
                dict(
                    bounds=np.column_stack((lo, hi)).tolist(),
                    medium={"density_kg_m3": 1000, "sound_speed_m_s": 1500, "viscosity_Pa_s": 0.002},
                )
            ),
        )
        h["lower"] = lo[None, :]
        h["upper"] = hi[None, :]
        h["valid"] = [True]
        h["harmonics"] = np.repeat(v[:, None], 8, axis=1)[None, :, :]
        h["max_R"] = np.full((1, 8), 1.1)
        h["min_R"] = np.full((1, 8), 0.9)
    return path


def test_adaptive_affine_complex_phase_and_invalid_cell(adaptive_fixture):
    m = AdaptiveResponseTableBubbleModel(adaptive_fixture)
    q = m.query(1.5e-6, 1.5e6, 1.5e4)
    np.testing.assert_allclose(q["harmonics_m3"], (3 + 1.5j) * 1e-18, atol=0, rtol=1e-14)
    np.testing.assert_allclose(
        m.harmonic_response(1.5e-6, 1.5e6, 1.5e4j)["harmonics_m3"],
        (3 + 1.5j) * 1e-18 * 1j ** np.arange(1, 9),
        atol=1e-32,
    )
    for r, f, p in [(1e-9, 1.5e6, 1.5e4), (1.5e-6, 9e6, 1.5e4), (1.5e-6, 1.5e6, np.nan)]:
        with pytest.raises(ValueError):
            m.query(r, f, p)
    with h5py.File(adaptive_fixture, "r+") as h:
        h["valid"][:] = False
    with pytest.raises(ValueError, match="INVALID_CELL"):
        AdaptiveResponseTableBubbleModel(adaptive_fixture).query(1.5e-6, 1.5e6, 1.5e4)
    with h5py.File(adaptive_fixture, "r+") as h:
        h.attrs["complete"] = False
    with pytest.raises(ValueError, match="Incomplete"):
        AdaptiveResponseTableBubbleModel(adaptive_fixture)


def test_gate_missing_scope_tamper_and_override(tmp_path):
    evidence = tmp_path / "evidence"
    evidence.write_text("complete")
    path = tmp_path / "gates.json"
    manager = ValidationGateManager(path)
    with pytest.raises(ValidationGateError):
        manager.check("run-lazarus", scope="paper")
    m = manager.check("run-lazarus", scope="paper", force_demo=True)
    assert m["UNVALIDATED_OVERRIDE"] and validation_watermark(m) == "DEMO / NOT VALIDATED"
    from mevo_mb.validation.gates import REQUIRES

    path.write_text(
        json.dumps(
            {
                "gates": {
                    k: dict(
                        status="PASS",
                        scopes=["paper"],
                        evidence=[
                            {
                                "path": str(evidence),
                                "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
                            }
                        ],
                    )
                    for k in REQUIRES["run-lazarus"]
                }
            }
        )
    )
    admitted = ValidationGateManager(path).check("run-lazarus", scope="paper")
    assert admitted["gate_status"] == "PASS" and "physics_status" not in admitted
    with pytest.raises(ValidationGateError):
        ValidationGateManager(path).check("run-lazarus", scope="sonovue")
    evidence.write_text("modified")
    with pytest.raises(ValidationGateError):
        ValidationGateManager(path).check("run-lazarus", scope="paper")


def test_measured_regime_rejects_unsupported_shell_frequency(tmp_path):
    rows = []
    for model in ["zhang", "doinikov"]:
        rows.append(
            dict(
                model=model,
                shell_model="uncoated",
                sound_speed_m_s=1482.0,
                radius1_m=1.5e-6,
                radius2_m=2e-6,
                separation_m=20e-6,
                frequency_Hz=300000.0,
                pressure_Pa=100.0,
                max_radial_fraction=0.001,
                relative_distance_change=0.0001,
                calculated_sign=1,
                classical_sign=1,
                coefficient_relative_error=0.01,
                hydrodynamic_relative_error=0.02,
            )
        )
    p = tmp_path / "map.json"
    p.write_text(json.dumps(dict(rows=rows)))
    m = PairPhysicsRegimeMap(p)
    args = [1.5e-6, 2e-6, 20e-6, 3e5, 100, "uncoated"]
    assert m.classify(*args) == "FAR_FIELD_SAFE"
    assert m.classify(*args[:-1], "sonovue") == "OUTSIDE_VALIDATED_RANGE"
    args[3] = 301000.0
    assert m.classify(*args) == "OUTSIDE_VALIDATED_RANGE"
    args[2] = 3e-6
    assert m.classify(*args) == "CONTACT"


def test_explicit_gpu_missing_and_oom_fallback(monkeypatch):
    from mevo_mb.bubbles.auto_backend import AutoPairBackend

    x = np.array([[0.0, 0, 0], [25e-6, 0, 0]])
    args = (x, np.ones(2) * 1e-6, np.ones(2) * 1e-18, 1e6, 1000.0)
    for reason in ["missing", "oom"]:
        model = AutoPairBackend(crossover=1)
        if reason == "missing":

            def unavailable():
                raise RuntimeError("GPU unavailable injected")

            monkeypatch.setattr(model, "_gpu", unavailable)
        else:

            class GPU:
                def forces(self, *a):
                    raise MemoryError("injected CUDA OOM")

            model.gpu = GPU()
        with pytest.warns(RuntimeWarning):
            result = model.forces(*args)
        np.testing.assert_array_equal(result, model.cpu.forces(*args))
        assert model.gpu_disabled and model.decisions[-1]["backend"] == "numpy"


@pytest.mark.parametrize(
    "kw", [{"pressure_Pa": np.nan}, {"radius_m": -1e-6}, {"elasticity_N_m": -1}, {"sigma_initial_N_m": 0.08}]
)
def test_bad_radial_inputs(kw):
    with pytest.raises(ValueError):
        replace(case(), **kw)


def test_concurrent_duplicate_cache_calls_are_serialized(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from mevo_mb.bubbles.apecss import ApecssBubbleResponseModel
    import time

    exe = tmp_path / "apecss"
    exe.write_bytes(b"fixture")
    cache = ResponseNodeCache(tmp_path / "cache", exe)
    calls = []

    def solve(self, c, out):
        calls.append(c)
        time.sleep(0.03)
        out.mkdir()
        np.savez(out / "cycles.npz", harmonics_cosine_m3=np.ones(8) * 1e-18, radius_m=np.array([c.radius_m]))
        meta = dict(status="COMPLETE", metrics={"valid": True}, returncode=0, elapsed_s=0.03)
        (out / "RESULT.json").write_text(json.dumps(meta))
        return meta, None, None

    monkeypatch.setattr(ApecssBubbleResponseModel, "solve", solve)
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(cache.get, [case()] * 4))
    assert len(calls) == 1 and sum(r["cache_hit"] for r in rows) == 3


def test_adaptive_builder_refinement_checkpoint_resume(tmp_path, monkeypatch):
    import mevo_mb.bubbles.adaptive_table as module

    class Pool:
        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def map(self, fn, tasks):
            return map(fn, tasks)

    calls = []

    def task(args):
        c = args[2]
        calls.append(c)
        h = np.zeros(8, complex)
        h[0] = 1e-18 * (1 + (c.pressure_Pa / 40000) ** 2) * (1 + 0.3j)
        return dict(
            status="COMPLETE",
            cache_hit=False,
            harmonics_real=h.real.tolist(),
            harmonics_imag=h.imag.tolist(),
            metrics=dict(
                valid=True,
                harmonic_energy_resolved=True,
                max_R_over_R0=1.05,
                min_R_over_R0=0.95,
                mean_dot_volume_squared_m6_s2=float(0.5 * abs(2 * np.pi * c.frequency_Hz * h[0]) ** 2),
            ),
        )

    monkeypatch.setattr(module, "ProcessPoolExecutor", Pool)
    monkeypatch.setattr(module, "cached_task", task)
    exe = tmp_path / "apecss"
    exe.write_bytes(b"fixture")
    args = (tmp_path / "tree", exe, Path("configs/stage2/sonovue_parameters.json"))
    first = module.AdaptiveResponseTableBuilder(*args, workers=1).build()
    assert first["coverage_fraction"] == pytest.approx(1) and first["leaves"] > 1
    n = len(calls)
    resumed = module.AdaptiveResponseTableBuilder(*args, workers=1).build()
    assert resumed["new_apecss_calls"] == 0 and len(calls) == n
    with pytest.raises(ValueError, match="mismatch"):
        module.AdaptiveResponseTableBuilder(*args, workers=1, profile="strict")
    tree = json.loads((args[0] / "CHECKPOINT.json").read_text())
    for leaf in tree["leaves"]:
        assert max(leaf["errors"].values()) <= 0.05


def test_adaptive_gpu_batch_matches_cpu(adaptive_fixture):
    cp = pytest.importorskip("cupy")
    try:
        count = cp.cuda.runtime.getDeviceCount()
    except Exception:
        pytest.skip("CUDA unavailable")
    if not count:
        pytest.skip("CUDA unavailable")
    cpu = AdaptiveResponseTableBubbleModel(adaptive_fixture)
    gpu = AdaptiveResponseTableBubbleModel(adaptive_fixture, backend="cupy")
    args = (np.linspace(1e-6, 2e-6, 5000), 1.5e6, 1.3e4)
    a = cpu.query(*args)
    b = gpu.query(*args)
    for key in ["harmonics_m3", "max_R_over_R0", "min_R_over_R0"]:
        np.testing.assert_allclose(cp.asnumpy(b[key]), a[key], rtol=3e-15, atol=1e-32)


def test_rejected_table_and_corrupt_radius(adaptive_fixture):
    with h5py.File(adaptive_fixture, "r+") as h:
        h["min_R"][:] = np.nan
    with pytest.raises(ValueError, match="Invalid adaptive cell"):
        AdaptiveResponseTableBubbleModel(adaptive_fixture)
    (adaptive_fixture.parent / "TABLE_REJECTION.json").write_text("{}")
    with pytest.raises(ValueError, match="explicitly rejected"):
        AdaptiveResponseTableBubbleModel(adaptive_fixture)


def test_closed_gate_precedes_acoustic_load(tmp_path):
    from mevo_mb.coupling.experiment import run_particles

    config = {
        "bubble_response": {"model": "adaptive_response_table"},
        "validation_gates": {"registry": str(tmp_path / "absent")},
    }
    with pytest.raises(ValidationGateError):
        run_particles(config, tmp_path / "nonexistent_basis", tmp_path)
    decision = json.loads((tmp_path / "GATE_DECISION.json").read_text())
    assert decision["gate_status"] == "NOT_RUN_GATE" and not decision["UNVALIDATED_OVERRIDE"]


def test_translation_map_rejects_unmeasured_shell(tmp_path):
    from mevo_mb.validation.regime_maps import TranslationValidityMap

    p = tmp_path / "translation.json"
    p.write_text(
        json.dumps(
            {
                "rows": [
                    dict(
                        radius_m=2e-6,
                        pressure_Pa=150e3,
                        frequency_Hz=1e6,
                        shell_model="solid15nm",
                        cycles=20,
                        status="FAST_MODEL_SAFE",
                    )
                ]
            }
        )
    )
    m = TranslationValidityMap(p)
    assert m.classify(2e-6, 150e3, 1e6, shell_model="solid15nm") == "FAST_MODEL_SAFE"
    assert m.classify(2e-6, 150e3, 1e6, shell_model="sonovue") == "OUTSIDE_VALIDATED_RANGE"
    assert m.classify(2e-6, 150e3, 1e6, shell_model="solid15nm", cycles=20000) == "OUTSIDE_VALIDATED_RANGE"


def test_monitor_records_unknown_pair_warning_without_promoting_physics():
    from mevo_mb.bubbles.validity import PhysicsValidityMonitor

    m = PhysicsValidityMonitor("adaptive_response_table", policy="strict")
    m.observe(
        np.array([[0.0, 0, 0], [25e-6, 0, 0]]),
        np.full(2, 1e-6),
        0.0,
        radial_fraction=np.zeros(2),
        frequency=1e6,
        sound_speed=1500.0,
        pressure_amplitude=np.full(2, 1e4),
        shell_model="sonovue",
    )
    result = m.summary()
    assert result["pair_regime"]["invalid_pair_fraction"] == 1
    assert result["validity_status"] == "PHYSICS_VALIDITY_WARNING"
    assert result["physics_status"] == "partially validated"

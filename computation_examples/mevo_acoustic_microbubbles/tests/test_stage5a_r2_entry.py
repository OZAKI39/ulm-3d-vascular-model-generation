from types import SimpleNamespace
import pytest
from scripts.stage5a_r2 import run_v4a as entry


@pytest.mark.parametrize("step,factor", [(5e-6, 0.25), (2.5e-6, 0.0625)])
def test_stricter_formal_entry_dispatches_unchanged_loop(monkeypatch, step, factor):
    monkeypatch.setattr(entry, "verify", lambda: dict(resources=dict(per_run_s=1800)))
    calls = []
    monkeypatch.setattr(entry, "frozen_run", lambda *a, **kw: calls.append((a, kw)))
    args = SimpleNamespace(
        case="B2",
        duration=0.001,
        max_dt=step,
        tol_factor=factor,
        pair="cupy",
        nonsmooth_recovery=True,
        nonlinear_strategy="quasi_newton",
        label="formal_v4a_test",
        budget=1800,
    )
    entry.execute(args)
    assert len(calls) == 1 and calls[0][1]["method"] == "S1" and calls[0][1]["components"]
    assert calls[0][1]["max_dt"] == step and calls[0][1]["tol_factor"] == factor
    assert calls[0][1]["label"] == "validated_v4a_test"


@pytest.mark.parametrize(
    "key,value",
    [
        ("pair", "numpy"),
        ("budget", 1801),
        ("tol_factor", 1.0),
        ("nonsmooth_recovery", False),
        ("nonlinear_strategy", "picard"),
    ],
)
def test_formal_entry_rejects_outside_protocol_before_run(monkeypatch, key, value):
    monkeypatch.setattr(entry, "verify", lambda: dict(resources=dict(per_run_s=1800)))

    def forbidden(*a, **kw):
        raise AssertionError("Run must not initialize")

    monkeypatch.setattr(entry, "frozen_run", forbidden)
    args = SimpleNamespace(
        case="B2",
        duration=0.001,
        max_dt=5e-6,
        tol_factor=0.25,
        pair="cupy",
        nonsmooth_recovery=True,
        nonlinear_strategy="quasi_newton",
        label="formal_v4a_test",
        budget=1800,
    )
    setattr(args, key, value)
    with pytest.raises(ValueError, match="CONFIG_ERROR"):
        entry.execute(args)

"""Post-campaign diagnostic wrappers; no frozen scientific kernel is edited.

Every nested timing subtracts child instrumentation scopes. Inclusive QP-native
timings remain a separate table, never added to the exclusive partition.
"""

import argparse
from collections import defaultdict
from contextlib import contextmanager
from functools import wraps
from time import perf_counter
from scripts.stage5a_r3.common import BASE, write, sha
from scripts.stage5a_r3 import run as driver
from mevo_mb.bubbles.motion_r3 import SemiImplicitR3
from mevo_mb.bubbles.motion_backend_r3 import MotionR3Optimized
from mevo_mb.bubbles.motion_qp_r3 import NativeQPR3
from mevo_mb.bubbles.motion_qp_r2 import NativeQP
from mevo_mb.bubbles.motion_sweep_r3 import WallSweepR3


class Accounting:
    def __init__(self):
        self.stack = []
        self.exclusive = defaultdict(float)
        self.calls = defaultdict(int)
        self.trials, self.stages, self.qps, self.native = [], [], [], []
        self.trial_id = -1
        self.stage = None

    @contextmanager
    def scope(self, label):
        record = [perf_counter(), 0.0]
        self.stack.append(record)
        self.calls[label] += 1
        try:
            yield
        finally:
            elapsed = perf_counter() - record[0]
            self.stack.pop()
            self.exclusive[label] += elapsed - record[1]
            if self.stack:
                self.stack[-1][1] += elapsed

    def bind(self, cls, name, label, table=None):
        original = getattr(cls, name)

        @wraps(original)
        def wrapped(obj, *args, **kwargs):
            begin = perf_counter()
            result, error = None, None
            try:
                with self.scope(label):
                    result = original(obj, *args, **kwargs)
                return result
            except Exception as exc:
                error = repr(exc)
                raise
            finally:
                if table is not None:
                    d = result[-1] if result is not None else None
                    table.append(
                        dict(
                            trial_id=self.trial_id,
                            stage=self.stage,
                            wall_s=perf_counter() - begin,
                            error=error,
                            diagnostics=d,
                        )
                    )

        setattr(cls, name, wrapped)
        return original


def main(window):
    from mevo_mb.bubbles.tracker import ParticleTracker
    from scripts.stage5a_r3.current import campaign as current_campaign

    # Must remain isolated from all formal trajectory timings.
    campaign = BASE / "final/CAMPAIGN_R3_V1/RESULT.json"
    from scripts.stage5a_r3.common import read

    if campaign.exists():
        assert current_campaign()["complete"], "Do not contaminate formal campaign timings"
    b3 = window.startswith("B3_")
    b3_fixture = None
    if b3:
        import h5py
        import numpy as np
        from pathlib import Path
        from mevo_mb.bubbles.constrained_motion import AcceptedCache

        ledger = read(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")["ledger"]
        item = next(
            r for r in ledger if r["case"] == "B3" and r["kind"] == "nominal" and r["duration_s"] == 0.005
        )
        path = Path(item["run"]) / "trajectories.h5"
        with h5py.File(path) as h:
            k = int(np.argmin(abs(h["time_s"][:] - 3.0015)))
            x, r, ids, mask, t = (
                h["positions_m"][k],
                h["radii_m"][:],
                h["particle_id"][:],
                h["active"][k],
                float(h["time_s"][k]),
            )
        b3_fixture = BASE / "fixtures/B3_costly_actual_3p0015.npz"
        if not b3_fixture.exists():
            meta = dict(
                source=str(path),
                source_sha256=sha(path),
                actual_index=k,
                cache="Explicit cold diagnostic start; every original active particle retained",
            )
            driver.checkpoint(b3_fixture, x, np.zeros_like(x), r, ids, mask, t, meta, AcceptedCache(), meta)
    audit = Accounting()
    patches = []
    for cls, name, label, table in [
        (ParticleTracker, "evaluate", "force_response_flow", None),
        (SemiImplicitR3, "wall", "wall_query", None),
        (SemiImplicitR3, "feature_constraints", "finite_feature_constraints", None),
        (MotionR3Optimized, "solve", "resistance_assembly", None),
        (MotionR3Optimized, "constrained_stage", "qp_partition_and_global_certificate", None),
        (NativeQPR3, "solve", "qp_adapter_and_fallback", audit.qps),
        (NativeQP, "candidate", "qp_native_candidate", audit.native),
        (WallSweepR3, "__call__", "continuous_geometry", None),
    ]:
        patches.append((cls, name, audit.bind(cls, name, label, table)))

    class MeasuredStep(SemiImplicitR3):
        def advance(self, x, r, ids, t, dt, stage):
            before = dict(audit.exclusive)
            before_qp_calls = self.calls
            start = perf_counter()
            audit.stage = stage
            result, error = None, None
            try:
                with audit.scope("nonlinear_and_substep_other"):
                    result = super().advance(x, r, ids, t, dt, stage)
                return result
            except Exception as exc:
                error = repr(exc)
                raise
            finally:
                audit.stages.append(
                    dict(
                        trial_id=audit.trial_id,
                        stage=stage,
                        method=self.method,
                        time_s=t,
                        dt_s=dt,
                        wall_s=perf_counter() - start,
                        error=error,
                        exclusive_s={k: v - before.get(k, 0) for k, v in audit.exclusive.items()},
                        diagnostics=result[-1] if result is not None else None,
                        attempted_nonlinear_QP_iterations=self.calls - before_qp_calls,
                    )
                )
                audit.stage = None

        def trial(self, x, r, ids, t, dt):
            audit.trial_id += 1
            begin = perf_counter()
            first = len(audit.stages)
            error, result = None, None
            try:
                with audit.scope("step_doubling_error_and_other"):
                    result = super().trial(x, r, ids, t, dt)
                return result
            except Exception as exc:
                error = repr(exc)
                raise
            finally:
                audit.trials.append(
                    dict(
                        trial_id=audit.trial_id,
                        method=self.method,
                        time_s=t,
                        dt_s=dt,
                        wall_s=perf_counter() - begin,
                        error=error,
                        accepted=error is None,
                        stage_indices=list(range(first, len(audit.stages))),
                        diagnostics=result[-1] if result is not None else None,
                    )
                )

    original_step = driver.SemiImplicit
    driver.SemiImplicit = MeasuredStep
    original_checkpoint = driver.checkpoint

    def measured_checkpoint(*args, **kwargs):
        with audit.scope("checkpoint_IO"):
            return original_checkpoint(*args, **kwargs)

    driver.checkpoint = measured_checkpoint
    fixture = (
        b3_fixture
        if b3
        else {
            "initial": None,
            "stable": BASE / "fixtures/B2_actual_3.00300000.npz",
            "S0": BASE / "diagnosis/reproduction_v1/nominal/terminal.npz",
            "event": BASE / "diagnosis/reproduction_v1/nominal/before_persistent_sweep.npz",
        }[window]
    )
    begin = perf_counter()
    try:
        with audit.scope("driver_initialization_IO_monitoring_other"):
            out, summary = driver.run(
                "B3" if b3 else "B2",
                15e-6 if b3 else 4e-7 if window == "S0" else 1e-5,
                method="S0" if window == "S0" else "S1",
                max_dt=50e-9 if window == "S0" else 5e-6,
                tol_factor=0.25,
                pair="cupy",
                budget=300,
                label="dynamic_audit_" + window,
                initial_dt=5e-6 if b3 else 1e-8,
                components=window == "B3_serial_costly",
                nonlinear_strategy="quasi_newton",
                nonsmooth_recovery=window != "S0",
                fixture=fixture,
                optimized=True,
            )
    finally:
        driver.SemiImplicit = original_step
        driver.checkpoint = original_checkpoint
        for cls, name, original in reversed(patches):
            setattr(cls, name, original)
    wall = perf_counter() - begin
    accepted = {r["trial_id"] for r in audit.trials if r["accepted"]}
    stage_totals = []
    for method in ["S0", "S1"]:
        for stage in ["full", "half1", "half2"]:
            for accept in [True, False]:
                rows = [
                    r
                    for r in audit.stages
                    if r["method"] == method and r["stage"] == stage and (r["trial_id"] in accepted) == accept
                ]
                stage_totals.append(
                    dict(
                        method=method,
                        stage=stage,
                        trial_accepted=accept,
                        count=len(rows),
                        wall_s=sum(r["wall_s"] for r in rows),
                        nonlinear_iterations=sum(r["attempted_nonlinear_QP_iterations"] for r in rows),
                        failed_calls=sum(r["error"] is not None for r in rows),
                    )
                )
    result = dict(
        window=window,
        run=str(out),
        summary=summary,
        diagnostic_wall_s=wall,
        exclusive_s=dict(audit.exclusive),
        function_calls=dict(audit.calls),
        accounted_s=sum(audit.exclusive.values()),
        stage_totals=stage_totals,
        trials=audit.trials,
        stages=audit.stages,
        qp_calls=audit.qps,
        native_candidates=audit.native,
        scientific_kernel_modified=False,
        instrumentation_source_sha256=sha(__file__),
        scope="Diagnostic instrumentation overhead included; never use for completion speedup. Exclusive scopes form one partition; stage and native timing tables overlap this partition.",
        IO_scope="Checkpoints separately timed. Remaining HDF5/log IO and initialization belong to driver_other; original summary provides inclusive IO and monitor overhead.",
    )
    assert abs(result["accounted_s"] - wall) < 0.01
    write(BASE / "diagnosis" / f"DYNAMIC_INSTRUMENTED_{window}.json", result)
    print(dict(window=window, status=summary["status"], wall_s=wall, stages=stage_totals), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--window",
        choices=["initial", "stable", "S0", "event", "B3_global_costly", "B3_serial_costly"],
        required=True,
    )
    main(p.parse_args().window)

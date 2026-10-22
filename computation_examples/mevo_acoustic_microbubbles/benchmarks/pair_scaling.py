"""Bounded synthetic kernel benchmark; NOT an aggregation physics experiment."""

from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time

import cupy as cp
import numpy as np

from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from mevo_mb.bubbles.contact import constrain_pair_velocities, project_contacts
from mevo_mb.bubbles.response_table import ResponseTableBubbleModel
from mevo_mb.coupling.geometry import Cylinder


def timed(fun, sync=lambda: None, repeats=7):
    fun()
    sync()
    values = []
    for _ in range(repeats):
        sync()
        s = time.perf_counter()
        fun()
        sync()
        values.append(time.perf_counter() - s)
    return float(np.median(values))


def main():
    root = Path.cwd()
    out = root / "results/stage2/performance" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True)
    rows = []
    rng = np.random.default_rng(20261001)
    for n in (32, 64, 128, 256, 512, 1000):
        side = int(np.ceil(n ** (1 / 3)))
        x = np.array(np.unravel_index(np.arange(n), (side, side, side))).T * 30e-6
        r = np.full(n, 1e-6)
        v = (rng.normal(size=(n, 8)) + 1j * rng.normal(size=(n, 8))) * 1e-19
        reference = ClassicalFarFieldBjerknes().forces(x, r, v, 1e6, 1000.0)
        for name, xp in (("numpy", np), ("cupy", cp)):
            xx, rr, vv = map(xp.asarray, (x, r, v))
            backend = ClassicalFarFieldBjerknes(name)
            sync = (lambda: cp.cuda.Stream.null.synchronize()) if name == "cupy" else (lambda: None)
            def fun():
                return backend.forces(xx, rr, vv, 1e6, 1000.0)
            force_time = timed(fun, sync)
            actual = fun()
            actual = actual.get() if name == "cupy" else actual
            error = float(np.linalg.norm(actual - reference) / np.linalg.norm(reference))
            def detect():
                return xp.sum(
                            xp.triu(
                                xp.linalg.norm(xx[:, None, :] - xx[None, :, :], axis=-1) < rr[:, None] + rr[None, :], 1
                            )
                        )
            detection = timed(detect, sync)
            # Production contact feasibility remains the tested CPU resolver.
            contact = timed(
                lambda: constrain_pair_velocities(x, r, np.ones_like(x) * 0.001, 1e-12, 100), repeats=3
            )
            cyl = Cylinder(0.01, 0.02)

            def step():
                f = fun()
                f = f.get() if name == "cupy" else f
                velocity = constrain_pair_velocities(x, r, f / (6 * np.pi * 0.001 * r[:, None]), 1e-12, 100)
                return project_contacts(x + 1e-6 * velocity, r, cyl, 1e-12, 100)

            step_time = timed(step, sync, repeats=3)
            rows.append(
                dict(
                    N=n,
                    number_of_pairs=n * (n - 1) // 2,
                    backend=name,
                    pair_force_s=force_time,
                    contact_detection_s=detection,
                    contact_constraint_cpu_s=contact,
                    synthetic_step_s=step_time,
                    relative_force_error=error,
                    gpu_memory_pool_reserved_bytes=int(cp.get_default_memory_pool().total_bytes())
                    if name == "cupy"
                    else 0,
                    contact_count=0,
                    case="spaced 3D lattice, 8 harmonics, synthetic step (not RK integration)",
                )
            )
            print(rows[-1], flush=True)
    # A bounded dense contact chain diagnoses resolver cost; no trajectory generated.
    dense = []
    for n in (32, 64, 128):
        x = np.zeros((n, 3))
        x[:, 0] = np.arange(n) * 2e-6
        r = np.full(n, 1e-6)
        vel = np.zeros_like(x)
        vel[:, 0] = np.linspace(0.001, -0.001, n)
        duration = timed(lambda: constrain_pair_velocities(x, r, vel, 1e-12, 100), repeats=3)
        dense.append(dict(N=n, contact_constraint_cpu_s=duration))
    state = json.loads((root / "results/stage2/single_bubble/STATE.json").read_text())
    table = root / state["table"]
    # Same exact validated node, varied global phase; then a validated cell midpoint.
    gpu = ResponseTableBubbleModel(table, backend="cupy")
    cpu = ResponseTableBubbleModel(table)
    cells = np.argwhere(cpu.cell_valid)
    idx = cells[0]
    q = [float((a[i] + a[i + 1]) / 2) for a, i in zip(cpu.axes, idx)]
    n = 10000
    r = np.full(n, q[0])
    pressure = np.full(n, q[2])
    a = cpu.evaluate_batch(r, q[1], pressure)
    b = gpu.evaluate_batch(cp.asarray(r), q[1], cp.asarray(pressure)).get()
    lookup = dict(
        N=n,
        query_SI=q,
        relative_error=float(np.linalg.norm(a - b) / np.linalg.norm(a)),
        numpy_s=timed(lambda: cpu.evaluate_batch(r, q[1], pressure)),
        cupy_s=timed(
            lambda: gpu.evaluate_batch(cp.asarray(r), q[1], cp.asarray(pressure)),
            lambda: cp.cuda.Stream.null.synchronize(),
        ),
        table=str(table.relative_to(root)),
    )
    summary = dict(
        status="PASS"
        if max(x["relative_force_error"] for x in rows) < 1e-11 and lookup["relative_error"] < 1e-11
        else "FAIL",
        rows=rows,
        dense_contact_probe=dense,
        response_lookup=lookup,
        seed=20261001,
        hardware="results/stage2/system/SERVER_DOCTOR.json",
        git_commit=subprocess.check_output(
            ["git", "-c", f"safe.directory={root}", "rev-parse", "HEAD"], text=True
        ).strip(),
        limitations=[
            "No unlimited N² allocation: max_bubbles=2000 by default",
            "GPU pool reserved bytes are allocator high-water, not precise device peak VRAM",
            "Contact projection/NNLS remains CPU; GPU only accelerates pair force and detection",
            "Synthetic array layout is a performance fixture, not a validated physical cloud",
        ],
    )
    (out / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out.parent / "STATE.json").write_text(
        json.dumps(dict(summary, run=str(out.relative_to(root))), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()

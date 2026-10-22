"""Whole implicit substeps on late development states, never holdout states."""

from concurrent.futures import ProcessPoolExecutor
import json
from time import perf_counter
import h5py
import numpy as np
from scripts.stage5a_r2.run import load
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write, sha
from mevo_mb.bubbles.motion_r2 import SemiImplicit
from mevo_mb.bubbles.motion_qp_r2 import NativeQP


def benchmark(task):
    case, path, fraction, layout = task
    tr, sweep, _, _, _ = load(case, "S1", "numpy", layout == "component")
    with h5py.File(path) as h:
        times = h["time_s"][:]
        target = times[0] + fraction * (times[-1] - times[0])
        index = int(np.argmin(abs(times - target)))
        mask = h["active"][index]
        x = h["positions_m"][index][mask]
        r = h["radii_m"][:][mask]
        ids = h["particle_id"][:][mask]
        t = float(times[index])
    if layout == "dense_global":
        tr.translation.qp_factory = lambda: NativeQP("prox_dense_bounded")
        tr.translation.qp = tr.translation.qp_factory()
    step = SemiImplicit(tr.translation, tr.cylinder, tr.evaluate, sweep, nonlinear_strategy="quasi_newton")
    rows = []
    out = BASE / "benchmarks/late_layout_v2"
    out.mkdir(exist_ok=True)
    for repeat in range(2):
        begin = perf_counter()
        try:
            y, u, cache, d = step.advance(x, r, ids, t, 5e-6, "late_development")
            status = "PASS"
            tr.translation.cache = cache
            error = None
            np.savez_compressed(
                out / f"{case}_{fraction:g}_{layout}_{repeat}.npz", x=x, y=y, u=u, ids=ids, radii_m=r
            )
        except Exception as exc:
            status = "FAIL"
            error = repr(exc)
            d = getattr(exc, "diagnostics", {})
        rows.append(
            dict(
                case=case,
                fraction=fraction,
                layout=layout,
                repeat=repeat,
                source=str(path),
                source_sha256=sha(path),
                frame=index,
                time_s=t,
                dt_s=5e-6,
                status=status,
                error=error,
                wall_s=perf_counter() - begin,
                diagnostics=d,
            )
        )
    return rows


if __name__ == "__main__":
    paths = {
        "B2": next(BASE.glob("runs/*B2_S1_0.001s_formal_0.001_nominal_repeat0/trajectories.h5")),
        "B3": next(BASE.glob("runs/*B3_S1_0.0001s_prefix_v4/trajectories.h5")),
    }
    tasks = [
        (c, p, f, layout)
        for c, p in paths.items()
        for f in [0.25, 0.5, 0.75, 0.99]
        for layout in ["sparse_global", "component", "dense_global"]
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=6) as pool:
        for records in pool.map(benchmark, tasks):
            rows.extend(records)
            print(
                json.dumps([{k: v for k, v in r.items() if k != "diagnostics"} for r in records]), flush=True
            )
    write(
        BASE / "benchmarks/late_layout_v2/RESULTS.json",
        dict(
            rows=rows,
            workers=6,
            threads=1,
            timing_scope="parallel diagnostic throughput; complete implicit substep including force/assembly/QP/geometry, not full trajectory latency",
            independent_holdout="fixtures/HOLDOUT_SELECTION_V2.json sealed before tuning; not opened here",
        ),
    )

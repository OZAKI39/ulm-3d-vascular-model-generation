"""Actual full/half1/half2 work at formal horizons and measured costly states."""

import subprocess
from pathlib import Path
from time import perf_counter, process_time
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, load, read, write, sha, identity
from mevo_mb.bubbles.motion_backend_r3 import MotionR3Optimized
from mevo_mb.bubbles.motion_r3 import SemiImplicitR3
from mevo_mb.bubbles.motion_sweep_r3 import WallSweepR3


def main():
    commands = subprocess.check_output(["ps", "-eo", "args"], text=True).splitlines()
    assert not any("-m scripts.stage5a_r3.run " in c for c in commands), "WAIT_FOR_ACTIVE_TRAJECTORY"
    root = BASE / "benchmarks/trial_layout"
    root.mkdir(exist_ok=False)
    campaign = read(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")["ledger"]
    started = perf_counter()
    rows = []
    comparisons = []
    for case, at in [("B2", 3.0015), ("B2", 3.0045), ("B3", 3.0015), ("B3", 3.0025), ("B3", 3.0049)]:
        entry = next(
            r for r in campaign if r["case"] == case and r["kind"] == "nominal" and r["duration_s"] == 0.005
        )
        path = Path(entry["run"]) / "trajectories.h5"
        with h5py.File(path) as h:
            j = int(np.searchsorted(h["time_s"][:], at, side="right") - 1)
            mask = h["active"][j]
            x0 = h["positions_m"][j][mask]
            r = h["radii_m"][:][mask]
            ids = h["particle_id"][:][mask]
            t0 = float(h["time_s"][j])
        name = f"{case}_{at:.7f}"
        np.savez_compressed(
            root / (name + "_input.npz"), positions_m=x0, radii_m=r, particle_id=ids, time_s=t0
        )
        trajectories = {}
        for layout, workers in [("global", 1), ("serial", 1), ("parallel2", 2), ("parallel4", 4)]:
            tr, sw, cfg, state, _ = load(case, "cupy", layout != "global", workers)
            back = MotionR3Optimized(
                cfg["medium"]["viscosity_Pa_s"], components=layout != "global", component_workers=workers
            )
            sweep = WallSweepR3(tr.cylinder, sw.path)
            step = SemiImplicitR3(
                back,
                tr.cylinder,
                tr.evaluate,
                sweep,
                nonlinear_strategy="quasi_newton",
                atol_R=0.0005,
                rtol=0.0025,
                gap_rtol=0.005,
            )
            x = x0.copy()
            t = t0
            positions = [x.copy()]
            for repeat in range(3):
                begin = perf_counter()
                cpu = process_time()
                before = dict(step.timings)
                calls = step.calls
                result = None
                error = None
                try:
                    states, d = step.trial(x, r, ids, t, 5e-6)
                    result = d
                except Exception as exc:
                    error = repr(exc)
                elapsed = perf_counter() - begin
                row = dict(
                    case=case,
                    time_s=t,
                    source_time_s=t0,
                    layout=layout,
                    workers=workers,
                    repeat=repeat,
                    macro_dt_s=5e-6,
                    status="PASS" if result is not None else "FAIL",
                    error=error,
                    wall_s=elapsed,
                    CPU_s=process_time() - cpu,
                    QP_calls=step.calls - calls,
                    exclusive_s={k: step.timings[k] - before[k] for k in step.timings},
                    diagnostics=result,
                    source=str(path),
                    source_sha256=sha(path),
                )
                rows.append(row)
                write(root / "PROGRESS.json", dict(rows=rows))
                if result is None:
                    break
                back.cache = states[-1][2]
                x = states[-1][0]
                t += 5e-6
                positions.extend([s[0].copy() for s in states])
            if back.executor:
                back.executor.shutdown()
            trajectories[layout] = np.array(positions)
            np.savez_compressed(
                root / (name + "_" + layout + ".npz"),
                positions_m=np.array(positions),
                particle_id=ids,
                radii_m=r,
            )
        for layout in ["serial", "parallel2", "parallel4"]:
            a, b = trajectories["global"], trajectories[layout]
            same = a.shape == b.shape
            err = float(np.linalg.norm(a - b, axis=2).max()) if same else None
            comparisons.append(
                dict(
                    case=case,
                    time_s=t0,
                    layout=layout,
                    max_position_error_m=err,
                    status="PASS" if same and err <= 1e-11 else "FAIL",
                )
            )
    result = dict(
        status="PASS" if all(r["status"] == "PASS" for r in rows + comparisons) else "FAIL",
        rows=rows,
        comparisons=comparisons,
        wall_s=perf_counter() - started,
        identity=identity(),
        scope="Complete adaptive trial work including discarded full step. Three changing macros with all active bubbles; real states from completed formal runs; cold workspace first macro included and distinguished.",
        physics_changed=False,
    )
    write(BASE / "benchmarks/LAYOUT_FULL_TRIAL.json", result)
    print(dict(status=result["status"], wall_s=result["wall_s"], rows=len(rows)), flush=True)


if __name__ == "__main__":
    main()

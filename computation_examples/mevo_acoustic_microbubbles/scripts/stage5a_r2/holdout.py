"""Open the six preregistered unused states only after source/test freeze."""

from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from time import perf_counter
import h5py
import numpy as np
from scipy import sparse
from scripts.stage5a_r2.freeze import verify
from scripts.stage5a_r2.seal import BASE
from scripts.stage5a_r2.run import load
from scripts.stage4l.dense_probe import sha, write
from mevo_mb.bubbles.motion_r2 import anchored_constraints, SemiImplicit
from mevo_mb.bubbles.constrained_motion import SubdivideInterval, MotionSolverFailure


def evaluate(item):
    verify()
    source = Path(item["source"])
    assert sha(source) == item["source_sha256"]
    with h5py.File(source) as h:
        k = item["index"]
        data = {
            key: h[key][k] if key in ["time_s", "positions_m", "active"] else h[key][:]
            for key in ["time_s", "positions_m", "active", "radii_m", "particle_id"]
        }
    raw = b"".join(np.asarray(v).tobytes() for v in data.values())
    assert hashlib.sha256(raw).hexdigest() == item["state_sha256"]
    tr, sw, _, _, _ = load(item["case"], "S1", "numpy", False)
    mask = data["active"]
    x = data["positions_m"][mask]
    r = data["radii_m"][mask]
    ids = data["particle_id"][mask]
    t = float(data["time_s"])
    out = BASE / "diagnosis/holdout_v2_component" / f"{item['case']}_{k}"
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(out / "SEALED_STATE.npz", **data)
    f = tr.evaluate(x, r, t)
    w = tr.cylinder.wall_gap_and_normal(x, r)
    a = tr.translation.solve(x, r, f["primary_force_N"] + f["secondary_force_N"], f["flow_velocity_m_s"], w)
    c = anchored_constraints(x, x, r, w, w, 1e-7, ids, 1e-11, a.velocity)
    sparse.save_npz(out / "A_SI.npz", a.matrix)
    sparse.save_npz(out / "G.npz", c.matrix)
    np.savez_compressed(
        out / "INPUT.npz",
        b=a.matrix @ a.velocity.ravel(),
        lower=c.lower,
        contact_ids=np.array([str(z) for z in c.ids]),
    )
    solutions = []
    audits = []
    for component, stage in [
        (False, "holdout_global_cold"),
        (False, "holdout_global_cold"),
        (True, "holdout_components"),
    ]:
        tr.translation.components = component
        u, cache, d = tr.translation.constrained_stage(a, c, ids, stage)
        tr.translation.cache = cache
        lam = np.array([cache.dual[key] for key in c.ids])
        reaction = c.matrix.T @ lam
        solutions.append((u, reaction))
        audits.append(d)
    tr.translation.components = True
    relative_u = [
        float(np.linalg.norm(u - solutions[0][0]) / max(np.linalg.norm(solutions[0][0]), 1e-300))
        for u, _ in solutions
    ]
    relative_reaction = [
        float(np.linalg.norm(v - solutions[0][1]) / max(np.linalg.norm(solutions[0][1]), 1e-300))
        for _, v in solutions
    ]
    step = SemiImplicit(tr.translation, tr.cylinder, tr.evaluate, sw, nonlinear_strategy="quasi_newton")
    elapsed = 0.0
    dt = 1e-8
    accepted = []
    rejected = []
    start = perf_counter()
    positions = [x.copy()]
    times = [t]
    try:
        while elapsed < 1e-6 - 1e-16:
            assert perf_counter() - start < 90 and len(rejected) < 100 and dt >= 1e-12, "HOLDOUT_BUDGET"
            dt = min(dt, 1e-6 - elapsed)
            try:
                states, d = step.trial(x, r, ids, t + elapsed, dt)
            except (SubdivideInterval, MotionSolverFailure) as exc:
                rejected.append(dict(dt=dt, elapsed=elapsed, diagnostics=exc.diagnostics))
                dt *= 0.5
                continue
            for j, (y, u, cache, sub) in enumerate(states):
                accepted.append(dict(time_s=t + elapsed + (j + 1) * dt / 2, dt_s=dt / 2, diagnostics=sub))
                positions.append(y.copy())
                times.append(t + elapsed + (j + 1) * dt / 2)
            x = states[-1][0]
            tr.translation.cache = states[-1][2]
            elapsed += dt
            dt *= min(2.0, max(0.5, 0.85 / max(d["error_ratio"], 1e-8) ** 0.5))
        status = "PASS" if max(relative_u) <= 1e-7 and max(relative_reaction) <= 1e-7 else "FAIL"
        error = None
    except Exception as exc:
        status = "FAIL"
        error = repr(exc)
    np.savez_compressed(out / "ACTUAL_MOTION.npz", time_s=times, positions_m=positions, ids=ids, radii_m=r)
    report = dict(
        selection=item,
        status=status,
        error=error,
        qp_audits=audits,
        relative_velocity_difference=relative_u,
        relative_reaction_difference=relative_reaction,
        completed_s=elapsed,
        accepted=accepted,
        rejected=rejected,
        wall_s=perf_counter() - start,
        source_protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_V2.json"),
        actual_motion_layout="component_serial",
        evaluator_sha256=sha(Path(__file__)),
    )
    write(out / "RESULT.json", report)
    return {
        k: report[k]
        for k in [
            "selection",
            "status",
            "error",
            "completed_s",
            "relative_velocity_difference",
            "relative_reaction_difference",
            "wall_s",
        ]
    }


if __name__ == "__main__":
    protocol = verify()
    selected = BASE / "fixtures/HOLDOUT_SELECTION_V2.json"
    assert sha(selected) == protocol["holdout_selection_sha256"]
    with ProcessPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(evaluate, json.loads(selected.read_text())))
    result = dict(
        status="PASS" if all(r["status"] == "PASS" for r in rows) else "FAIL",
        rows=rows,
        workers=6,
        BLAS_threads=1,
        pair_backend="numpy",
        timing="parallel diagnostic throughput, not formal latency",
    )
    write(BASE / "final/HOLDOUT_V2_COMPONENT.json", result)
    print(json.dumps(result), flush=True)

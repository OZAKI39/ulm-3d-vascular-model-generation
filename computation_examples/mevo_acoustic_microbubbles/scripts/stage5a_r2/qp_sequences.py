"""Real changing-state QP replay; 6 independent CPU workers, BLAS=1."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter
import json
import h5py
import numpy as np
from scipy import sparse
from scripts.stage4l.dense_probe import load_case, write, sha
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.bubbles.motion_r2 import anchored_constraints
from mevo_mb.bubbles.motion_qp_r2 import NativeQP

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver_r2"


def generate():
    out = BASE / "fixtures/development_sequence"
    out.mkdir(exist_ok=False)
    selection = json.loads(
        (ROOT / "results/stage5/motion_solver/audit_v2/BASELINE_SELECTION.json").read_text()
    )
    records = []
    for case in ["B2", "B3"]:
        cfg, field, flow, _, inputs = load_case(Path("/workspace/mevo_acoustic_microbubbles"), case)
        cfg["particles"].update(
            translation_backend="constrained_lubrication_r2", pair_backend="numpy", response_backend="numpy"
        )
        cfg["particles"]["constrained_solver"] = dict(components=False, solver="prox")
        tracker = ParticleTracker(cfg, field, flow)
        source = Path(next(r["run"] for r in selection if r["key"] == case + "_5ms")) / "trajectories.h5"
        with h5py.File(source) as h:
            for frame in [0, 1, 2, 3, 4]:
                mask = h["active"][frame]
                x = h["positions_m"][frame][mask]
                r = h["radii_m"][:][mask]
                ids = h["particle_id"][:][mask]
                t = float(h["time_s"][frame])
                f = tracker.evaluate(x, r, t)
                w = tracker.cylinder.wall_gap_and_normal(x, r)
                m = tracker.translation
                a = m.solve(x, r, f["primary_force_N"] + f["secondary_force_N"], f["flow_velocity_m_s"], w)
                c = anchored_constraints(x, x, r, w, w, 1e-7, ids, 1e-11, a.velocity)
                vs = float(np.max(abs(a.velocity)))
                an = a.matrix / a.scale
                bn = an @ (a.velocity.ravel() / vs)
                prefix = out / f"{case}_{frame:02d}"
                sparse.save_npz(str(prefix) + "_A.npz", an)
                sparse.save_npz(str(prefix) + "_G.npz", c.matrix)
                np.savez_compressed(str(prefix) + "_data.npz", b=bn, lower=c.lower / vs, ids=ids, x=x, r=r)
                records.append(
                    dict(
                        case=case,
                        frame=frame,
                        time_s=t,
                        source=str(source),
                        source_sha256=sha(source),
                        velocity_scale_m_s=vs,
                        resistance_scale_N_s_m=a.scale,
                        prefix=str(prefix),
                        contact_ids=[list(k) for k in c.ids],
                    )
                )
    write(out / "MANIFEST.json", records)


def replay(task):
    case, solver, mode = task
    manifest = json.loads((BASE / "fixtures/development_sequence/MANIFEST.json").read_text())
    seq = [r for r in manifest if r["case"] == case]
    rows = []
    q = NativeQP(solver, 1e-12 if solver != "clarabel" else 1e-10)
    prior_u = prior_l = prior_keys = None
    for index in range(5):
        item = seq[index if mode == "adjacent" else 0]
        p = item["prefix"]
        a, g = sparse.load_npz(p + "_A.npz"), sparse.load_npz(p + "_G.npz")
        d = dict(np.load(p + "_data.npz"))
        b = d["b"].copy()
        lo = d["lower"]
        keys = [tuple(k) for k in item["contact_ids"]]
        if mode == "perturbed":
            b *= 1 + index * 1e-5
        warm_u = prior_u if mode != "cold" else None
        warm_l = (
            np.array([dict(zip(prior_keys, prior_l)).get(k, 0.0) for k in keys])
            if prior_keys is not None and mode != "cold"
            else None
        )
        # Warm velocities/duals transformed through their physical scales when
        # real adjacent states use a different frozen normalization.
        if mode == "adjacent" and prior_u is not None:
            warm_u = warm_u * seq[index - 1]["velocity_scale_m_s"] / item["velocity_scale_m_s"]
            warm_l = (
                warm_l
                * (seq[index - 1]["velocity_scale_m_s"] * seq[index - 1]["resistance_scale_N_s_m"])
                / (item["velocity_scale_m_s"] * item["resistance_scale_N_s_m"])
            )
        if mode == "cold":
            q = NativeQP(solver, q.eps)
        start = perf_counter()
        u, l, audit = q.candidate(a, b, g, lo, warm_u, warm_l)
        row = dict(
            case=case,
            solver=solver,
            mode=mode,
            index=index,
            frame=item["frame"],
            time_s=item["time_s"],
            wall_s=perf_counter() - start,
            diagnostics=audit,
        )
        folder = BASE / "benchmarks/changing_qp"
        folder.mkdir(exist_ok=True)
        target = folder / f"{case}_{solver}_{mode}_{index}.npz"
        np.savez_compressed(target, u=u, lambda_=l, reaction=g.T @ l)
        row["candidate_file"] = str(target)
        rows.append(row)
        prior_u, prior_l, prior_keys = u, l, keys
    return rows


if __name__ == "__main__":
    generate()
    tasks = [
        (c, s, m)
        for c in ["B2", "B3"]
        for s in ["prox", "osqp", "clarabel"]
        for m in ["cold", "same", "adjacent", "perturbed"]
    ]
    results = []
    with ProcessPoolExecutor(max_workers=6) as pool:
        for rows in pool.map(replay, tasks):
            results.extend(rows)
            print(
                rows[0]["case"],
                rows[0]["solver"],
                rows[0]["mode"],
                max(r["diagnostics"]["kkt"] for r in rows),
                flush=True,
            )
    write(
        BASE / "benchmarks/changing_qp/RESULTS.json",
        dict(
            rows=results,
            resource_strategy="6 independent CPU diagnostic workers; BLAS/OpenMP=1; parallel throughput, not isolated latency",
        ),
    )

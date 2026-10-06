"""Recorded failed holdouts become explicit diagnostic inputs; no tuning in this replay."""

from concurrent.futures import ProcessPoolExecutor
import json
import numpy as np
from scipy.sparse.linalg import norm as sparse_norm
from scripts.stage5a_r2.run import load
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write, sha
from mevo_mb.bubbles.motion_r2 import SemiImplicit


def replay(task):
    index, dt, method = task
    folder = BASE / "diagnosis/holdout_v2_component" / f"B2_{index}"
    state = np.load(folder / "ACTUAL_MOTION.npz")
    x = state["positions_m"][-1]
    r = state["radii_m"]
    ids = state["ids"]
    t = float(state["time_s"][-1])
    tr, sw, _, _, _ = load("B2", method, "numpy", True)
    model = tr.translation
    original_solve = model.solve
    original_constrained = model.constrained_stage
    snapshots = []
    qp = []
    matrices = []
    arrays = {}

    def traced(position, radius, force, flow, wall):
        result = original_solve(position, radius, force, flow, wall)
        i = len(snapshots)
        normal_change = (
            np.linalg.norm(wall[1] - arrays[f"wall_normal_{i - 1}"], axis=1) if i else np.zeros(len(r))
        )
        worst = int(np.argmax(normal_change))
        snapshots.append(
            dict(
                index=i,
                normal_change_max=float(normal_change.max()),
                normal_change_particle_ID=int(ids[worst]),
                gap_at_largest_normal_change_m=float(wall[0][worst]),
                wall_coefficient_at_change=float(result.wall_coefficients_N_s_m[worst]),
                A_relative_change=float(
                    sparse_norm(result.matrix - matrices[-1]) / sparse_norm(result.matrix)
                )
                if i
                else 0.0,
                largest_position_change_m=float(np.linalg.norm(position - arrays[f"x_{i - 1}"], axis=1).max())
                if i
                else 0.0,
                pair_count=len(result.pairs),
            )
        )
        arrays[f"x_{i}"] = position.copy()
        arrays[f"wall_normal_{i}"] = wall[1].copy()
        arrays[f"wall_gap_{i}"] = wall[0].copy()
        arrays[f"pair_ids_{i}"] = ids[result.pairs]
        arrays[f"pair_normals_{i}"] = result.normals.copy()
        arrays[f"pair_gaps_{i}"] = result.gaps.copy()
        arrays[f"pair_coeff_{i}"] = (
            result.coefficients_N_s_m if hasattr(result, "coefficients_N_s_m") else np.array([])
        )
        matrices.append(result.matrix.copy())
        return result

    def traced_qp(*args, **kw):
        answer = original_constrained(*args, **kw)
        u, cache, d = answer
        qp.append(
            dict(
                kkt=d["kkt"],
                relative_force_residual=d["physical_force_residual_relative"],
                u_max=float(abs(u).max()),
            )
        )
        arrays[f"u_{len(qp) - 1}"] = u.copy()
        return answer

    model.solve = traced
    model.constrained_stage = traced_qp
    step = SemiImplicit(model, tr.cylinder, tr.evaluate, sw, method=method, nonlinear_strategy="quasi_newton")
    try:
        _, _, _, d = step.advance(x, r, ids, t, dt, "diagnostic")
        status = "PASS"
        error = None
    except Exception as exc:
        status = "FAIL"
        error = repr(exc)
        d = getattr(exc, "diagnostics", {})
    out = BASE / "diagnosis/nonlinear_v4"
    out.mkdir(exist_ok=True)
    label = f"B2_{index}_{method}_dt{dt:g}"
    np.savez_compressed(out / (label + ".npz"), initial_x=x, radii=r, ids=ids, time_s=t, **arrays)
    record = dict(
        status=status,
        error=error,
        source=str(folder / "ACTUAL_MOTION.npz"),
        source_sha256=sha(folder / "ACTUAL_MOTION.npz"),
        index=index,
        method=method,
        dt_s=dt,
        time_s=t,
        diagnostics=d,
        iterations=snapshots,
        qp=qp,
    )
    write(out / (label + ".json"), record)
    return dict(
        index=index,
        method=method,
        dt_s=dt,
        status=status,
        error=error,
        max_normal_change=max(r["normal_change_max"] for r in snapshots),
        last_records=snapshots[-4:],
    )


if __name__ == "__main__":
    tasks = [
        (i, d, m)
        for i in [17384, 15472]
        for d, m in [(1e-8, "S1"), (1e-10, "S1"), (1e-12, "S1"), (1e-10, "S0")]
    ]
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(replay, tasks))
    write(BASE / "diagnosis/nonlinear_v4/SUMMARY.json", rows)
    print(json.dumps(rows), flush=True)

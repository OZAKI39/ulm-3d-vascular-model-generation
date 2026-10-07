"""Reconstruct complete runs' final original equations from atomic checkpoints.

Run only after isolated timing jobs finish. NumPy force evaluation is an additional
CPU/GPU cross-check; original recorded online certificates are kept separately.
"""

from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.spatial.distance import cdist
from scripts.stage5a_r2.run import load, restore
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write, sha
from mevo_mb.bubbles.motion_r2 import anchored_constraints
from mevo_mb.bubbles.motion_qp_r2 import certificate


def audit(path):
    path = Path(path)
    summary = json.loads((path / "SUMMARY.json").read_text())
    state, meta, cache = restore(path / "CHECKPOINT.npz")
    with h5py.File(path / "trajectories.h5") as h:
        t0, t1 = map(float, h["time_s"][-2:])
        mask = h["active"][-1]
        initial = h["positions_m"][0][mask]
        x0 = h["positions_m"][-2][mask]
        x1 = h["positions_m"][-1][mask]
    ledger = [json.loads(s) for s in (path / "ACCEPTED.jsonl").read_text().splitlines()]
    last = ledger[-1]
    method = last.get("substep_method", "S1")
    dt = last["dt_s"]
    tr, sweep, _, _, _ = load(summary["case"], "S1", "numpy", True)
    r = state["radii_m"][mask]
    ids = state["particle_id"][mask]
    u = state["velocity_m_s"][mask]
    assert np.array_equal(x1, state["positions_m"][mask])
    f = tr.evaluate(x0, r, t0)
    force = f["primary_force_N"] + f["secondary_force_N"]
    flow = f["flow_velocity_m_s"]
    wall0 = tr.cylinder.wall_gap_and_normal(x0, r)
    a0 = tr.translation.solve(x0, r, force, flow, wall0)
    constraints = anchored_constraints(x0, x0, r, wall0, wall0, dt, ids, 1e-11, a0.velocity)
    position = x1 if method == "S1" else x0
    matrix = tr.translation.solve(position, r, force, flow, tr.cylinder.wall_gap_and_normal(position, r))
    vs = max(float(np.max(abs(matrix.velocity))), 1e-30)
    lam = np.array([cache.dual.get(k, 0.0) for k in constraints.ids]) / (matrix.scale * vs)
    an = matrix.matrix / matrix.scale
    bn = an @ (matrix.velocity.ravel() / vs)
    proof = certificate(
        an,
        bn,
        constraints.matrix,
        constraints.lower / vs,
        u.ravel() / vs,
        lam,
        velocity_scale=vs,
        force_scale=matrix.scale * vs,
        dt=dt,
    )
    safety = sweep(x0, x1, r, 1e-11)
    gaps = {}
    for name, points in [("initial", initial), ("final", x1)]:
        pg = cdist(points, points) - r[:, None] - r[None, :]
        np.fill_diagonal(pg, np.inf)
        wg = tr.cylinder.wall_gap_and_normal(points, r)[0]
        gaps[name] = dict(
            nearest_pair_signed_gap_quantiles_m=np.quantile(pg.min(0), [0, 0.25, 0.5, 0.75, 1]).tolist(),
            wall_signed_gap_quantiles_m=np.quantile(wg, [0, 0.25, 0.5, 0.75, 1]).tolist(),
            scope="All active original IDs, raw signed gaps; negative values preserved.",
        )
    return dict(
        run=str(path),
        method=method,
        checkpoint_sha256=sha(path / "CHECKPOINT.npz"),
        time_s=t1,
        dt_from_ledger_s=dt,
        timestamp_difference_s=t1 - t0,
        kinematic_roundoff_m=float(np.max(abs((x0 + dt * u) - x1))),
        certificate=proof,
        safety=safety,
        recorded_online_original_kkt=last["nonlinear_kkt"] if method == "S1" else last["kkt"],
        gap_distributions=gaps,
        pass_gate=bool(proof["pass_gate"] and safety["pass_gate"]),
        scope="Final actual substep only; independently rebuilt original matrix and all constraint rows. S0 audits start matrix; S1 audits endpoint matrix.",
    )


if __name__ == "__main__":
    campaign = json.loads(sorted((BASE / "final").glob("CAMPAIGN_V*/RESULT.json"))[-1].read_text())
    assert campaign["complete"]
    paths = [
        e["run"]
        for e in campaign["ledger"]
        if e.get("summary") and e["kind"] != "reference" and e["summary"]["status"] == "COMPLETE"
    ]
    with ProcessPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(audit, paths))
    write(
        BASE / "final/INDEPENDENT_ENDPOINT_AUDIT.json",
        dict(
            status="PASS" if all(r["pass_gate"] for r in rows) and rows else "FAIL",
            rows=rows,
            caveat="This supplements all online per-step certificates; it is not a recomputation of every accepted step.",
        ),
    )
    print([(Path(r["run"]).name, r["certificate"]["kkt"], r["pass_gate"]) for r in rows], flush=True)

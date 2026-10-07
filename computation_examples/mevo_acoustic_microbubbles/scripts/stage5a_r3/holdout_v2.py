"""New source-bound full-system solver holdouts; not trajectory truth labels."""

from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, load, read, write, sha
from scripts.stage5a_r3.freeze_v2 import verify
from scripts.stage5a_r3.run import run, checkpoint
from mevo_mb.bubbles.motion_r3 import MotionR3, SemiImplicitR3
from mevo_mb.bubbles.motion_backend_r3 import MotionR3Optimized
from mevo_mb.bubbles.motion_sweep_r3 import WallSweepR3
from mevo_mb.bubbles.constrained_motion import AcceptedCache


def main():
    protocol = verify()
    selected = read(BASE / "fixtures/HOLDOUT_R3_V2.json")
    rows = []
    for number, item in enumerate(selected):
        source = Path(item["source"])
        assert sha(source) == item["source_sha256"]
        with h5py.File(source) as h:
            k = item["index"]
            x, r, ids, mask, t = (
                h["positions_m"][k],
                h["radii_m"][:],
                h["particle_id"][:],
                h["active"][k],
                float(h["time_s"][k]),
            )
        fixture = BASE / "fixtures" / f"HOLDOUT_R3_V2_{number}.npz"
        checkpoint(
            fixture,
            x,
            np.zeros_like(x),
            r,
            ids,
            mask,
            t,
            dict(scope="cold solver holdout", source=item),
            AcceptedCache(),
            item,
        )
        tracker, sweep0, cfg, _, _ = load(item["case"], "cupy", False)
        back = MotionR3Optimized(cfg["medium"]["viscosity_Pa_s"], components=False)
        serial = MotionR3(cfg["medium"]["viscosity_Pa_s"], components=True)
        sweep = WallSweepR3(tracker.cylinder, sweep0.path)
        step = SemiImplicitR3(
            back, tracker.cylinder, tracker.evaluate, sweep, nonlinear_strategy="quasi_newton"
        )
        xx, rr, ii = x[mask], r[mask], ids[mask]
        f = tracker.evaluate(xx, rr, t)
        wall = tracker.cylinder.wall_gap_and_normal(xx, rr)
        a = back.solve(xx, rr, f["primary_force_N"] + f["secondary_force_N"], f["flow_velocity_m_s"], wall)
        c = step.feature_constraints(xx, rr, wall, 50e-9, ii, a.velocity)
        u, cache, d = back.constrained_stage(a, c, ii, "holdout")
        v, other, e = serial.constrained_stage(a, c, ii, "holdout")
        reaction = c.matrix.T @ np.array([cache.dual[key] for key in c.ids])
        reaction_old = c.matrix.T @ np.array([other.dual[key] for key in c.ids])
        dv = float(np.linalg.norm(u - v) / max(np.linalg.norm(v), 1e-30))
        dr = float(np.linalg.norm(reaction - reaction_old) / max(np.linalg.norm(reaction_old), 1e-30))
        layout = protocol["method"]["layouts"][item["case"]]
        out, summary = run(
            item["case"],
            500e-9,
            max_dt=5e-6,
            tol_factor=0.25,
            pair="cupy",
            budget=300,
            label=f"holdout_R3_V2_{number}",
            components=layout["components"],
            workers=layout["component_workers"],
            nonlinear_strategy="quasi_newton",
            nonsmooth_recovery=True,
            fixture=fixture,
            optimized=True,
        )
        numeric = (
            summary["status"] == "COMPLETE"
            and summary["accepted_kkt_max"] <= 1e-9
            and summary["accepted_linear_residual_max"] <= 1e-11
            and summary["projection_max_m"] == 0
            and summary["accepted_pair_swept_gap_min_m"] >= -1e-11
            and summary["accepted_wall_swept_gap_bound_min_m"] >= -1e-11
        )
        row = dict(
            selection=item,
            run=str(out),
            summary=summary,
            common_qp=dict(
                global_kkt=d["kkt"], serial_kkt=e["kkt"], velocity_relative=dv, reaction_relative=dr
            ),
            status="PASS" if numeric and dv <= 1e-7 and dr <= 1e-7 else "FAIL",
        )
        rows.append(row)
        write(BASE / "final/HOLDOUT_R3_V2.json", dict(status="RUNNING", rows=rows))
    status = "PASS" if all(r["status"] == "PASS" for r in rows) else "FAIL"
    write(
        BASE / "final/HOLDOUT_R3_V2.json",
        dict(status=status, rows=rows, protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_R3_V2.json")),
    )
    print(dict(status=status, count=len(rows)), flush=True)


if __name__ == "__main__":
    main()

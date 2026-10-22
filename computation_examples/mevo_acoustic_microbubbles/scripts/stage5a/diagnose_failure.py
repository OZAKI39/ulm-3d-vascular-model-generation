"""Bounded single-interval B2 diagnostic; no trajectory continuation or repair."""

from pathlib import Path
from datetime import datetime, timezone
import json
import traceback
import numpy as np
import h5py
from scipy.sparse import save_npz
from scripts.stage4l.dense_probe import load_case, write, sha
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.bubbles.constrained_motion import PredictorCorrector

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


def main():
    candidates = []
    for p in (BASE / "runs").glob("*/SUMMARY.json"):
        s = json.load(open(p))
        if s["case"] == "B2" and s["macro_dt_s"] == 20e-6 and s["status"] == "FAILED":
            candidates.append(p)
    p = sorted(candidates)[0]
    folder = p.parent
    identity = json.load(open(folder / "IDENTITY.json"))
    out = BASE / "diagnostics" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_failure_replay")
    out.mkdir(parents=True, exist_ok=False)
    config, field, flow, _, _ = load_case(Path("/workspace/mevo_acoustic_microbubbles"), "B2")
    config = identity["config"]
    tracker = ParticleTracker(config, field, flow)
    step = PredictorCorrector(
        tracker.translation,
        tracker.cylinder,
        tracker.evaluate,
        contact_tolerance_m=config["particles"]["contact_tolerance_m"],
    )
    with h5py.File(folder / "trajectories.h5") as h:
        mask = h["active"][-1]
        r = h["radii_m"][:][mask]
        ids = h["particle_id"][:][mask]
        x = h["positions_m"][-1][mask]
        previous = h["positions_m"][-2][mask]
        t = float(h["time_s"][-1])
        tp = float(h["time_s"][-2])
    # Reconstruct one accepted cache from the preceding original state. This is
    # an independent replay, not a bitwise saved-cache restart.
    record = dict(
        source_run=str(folder),
        source_summary_sha256=sha(p),
        source_trajectory_sha256=sha(folder / "trajectories.h5"),
        status="RUNNING",
        scope="Single-interval diagnostic. No new accepted trajectory; no parameter tuning.",
        warm_cache="reconstructed from one previous accepted interval, not saved exact cache",
        last_accepted_time_s=t,
    )
    try:
        _, _, cache, _ = step.trial(previous, r, ids, tp, t - tp)
        tracker.translation.cache = cache
        # ERROR.txt final recursive path LLLRRLRRLL: 10 binary splits of the
        # original 20 microsecond macro interval. Same physical endpoint formula.
        begin = 3.0
        end = 3.0 + 20e-6
        for direction in "LLLRRLRRLL":
            mid = begin + (end - begin) / 2
            if direction == "L":
                end = mid
            else:
                begin = mid
        if abs(begin - t) > 1e-15:
            raise ValueError("Failed interval reconstruction mismatch")
        record.update(attempt_dt_s=end - t, attempt_end_s=end)
        step.trial(x, r, ids, t, end - t)
        record["status"] = "FAILURE_NOT_REPRODUCED_WITH_RECONSTRUCTED_CACHE"
    except Exception as exc:
        record.update(
            status="REPLAY_REJECTED",
            exception=type(exc).__name__,
            diagnostics=getattr(exc, "diagnostics", None),
        )
        (out / "ERROR.txt").write_text(traceback.format_exc())
        trace = exc.__traceback__
        qp = None
        while trace:
            f = trace.tb_frame
            if f.f_code.co_name == "solve" and "saddle" in f.f_locals:
                qp = f.f_locals
            trace = trace.tb_next
        if qp:
            a, g = qp["a"], qp["g"]
            save_npz(out / "A.npz", a)
            save_npz(out / "G.npz", g)
            data = {
                key: qp[key]
                for key in [
                    "b",
                    "lower",
                    "u0",
                    "v",
                    "mult",
                    "working",
                    "direction",
                    "slack",
                    "warm_primal",
                    "warm_dual",
                    "feasible",
                ]
                if key in qp and qp[key] is not None
            }
            np.savez_compressed(out / "FAILED_QP_STATE.npz", **data)
            # Dense SVD is explicitly diagnostic-only on the small failed
            # component (114 columns), never a production contact path.
            singular = np.linalg.svd(g.toarray(), compute_uv=False)
            saddle = qp["saddle"]
            sv = np.linalg.svd(saddle.toarray(), compute_uv=False)
            dense = g.toarray()
            length = np.linalg.norm(dense, axis=1)
            cos = (dense / length[:, None]) @ (dense / length[:, None]).T
            np.fill_diagonal(cos, 0)
            maximum = np.max(cos, axis=1)
            record["matrix_audit"] = dict(
                primal_dimension=a.shape[0],
                constraint_count=g.shape[0],
                G_nnz=g.nnz,
                G_singular_values=singular.tolist(),
                active_saddle_shape=list(saddle.shape),
                active_saddle_smallest_singular=float(sv[-1]),
                active_saddle_condition=float(sv[0] / sv[-1]) if sv[-1] > 0 else None,
                rows_with_parallel_partner_cos_above_1_minus_1e_8=int(np.count_nonzero(maximum > 1 - 1e-8)),
                primal_fallback_iterations=int(qp.get("count", -1)) + 1,
                last_direction_norm=float(np.linalg.norm(qp["direction"])),
                reported_kkt_caveat="Failure-path summary computes stationarity with zero multipliers after fallback exhaustion; it is not an accepted-state KKT residual. Actual final working-set multiplier audit below.",
            )
            working = qp["working"]
            mult = qp["mult"]
            v = qp["v"]
            c = qp["c"]
            bounds = qp["bound"]
            if len(working) == len(mult):
                reaction = c[working].T @ mult
                stat = np.linalg.norm(a @ v - qp["b"] - reaction) / max(
                    1.0, np.linalg.norm(qp["b"]), np.linalg.norm(reaction)
                )
                record["matrix_audit"]["actual_last_working_stationarity"] = float(stat)
                record["matrix_audit"]["actual_last_working_min_multiplier"] = float(np.min(mult, initial=0))
                record["matrix_audit"]["actual_last_working_max_active_slack"] = float(
                    np.max(np.abs((c @ v - bounds)[working]), initial=0)
                )
    write(out / "DIAGNOSTIC.json", record)
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()

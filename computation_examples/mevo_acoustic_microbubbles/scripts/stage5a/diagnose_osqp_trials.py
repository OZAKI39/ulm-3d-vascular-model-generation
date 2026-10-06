"""Replay each rejected OSQP interval once, retaining its sparse QP inputs."""

from datetime import datetime, timezone
from pathlib import Path
import json
import traceback
import numpy as np
from scipy import sparse
from scripts.stage4l.dense_probe import load_case, sha, write
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.bubbles.constrained_motion import PredictorCorrector

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


def main():
    for folder in sorted((BASE / "runs").glob("*_osqp")):
        summary = json.loads((folder / "SUMMARY.json").read_text())
        if summary["accepted_steps"] or summary["status"] != "FAILED":
            continue
        out = (
            BASE
            / "diagnostics"
            / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_") + summary["case"] + "_osqp")
        )
        out.mkdir(parents=True, exist_ok=False)
        config = json.loads((folder / "CONFIG.json").read_text())
        _, field, flow, _, _ = load_case(Path("/workspace/mevo_acoustic_microbubbles"), summary["case"])
        tracker = ParticleTracker(config, field, flow)
        step = PredictorCorrector(
            tracker.translation,
            tracker.cylinder,
            tracker.evaluate,
            contact_tolerance_m=config["particles"]["contact_tolerance_m"],
        )
        with np.load(folder / "LAST_ACCEPTED_STATE.npz") as h:
            mask = h["active"]
            x, r, ids = h["positions_m"][mask], h["radii_m"][mask], h["particle_id"][mask]
            t = float(h["time_s"])
        rejected = [json.loads(line) for line in (folder / "REJECTED.jsonl").open()]
        last = rejected[-1]
        assert last["time_s"] == t
        end = t + (last["end_s"] - t) / 2
        record = dict(
            source_run=str(folder),
            source_summary_sha256=sha(folder / "SUMMARY.json"),
            case=summary["case"],
            time_s=t,
            end_s=end,
            dt_s=end - t,
            scope="One rejected interval; no accepted trajectory continuation",
            source_sha256={
                str(p.relative_to(ROOT)): sha(p)
                for p in [
                    Path(__file__),
                    ROOT / "src/mevo_mb/bubbles/osqp_motion.py",
                    ROOT / "src/mevo_mb/bubbles/constrained_motion.py",
                ]
            },
        )
        try:
            step.trial(x, r, ids, t, end - t)
            record["status"] = "FAILURE_NOT_REPRODUCED"
        except Exception as exc:
            record.update(
                status="REPLAY_REJECTED",
                exception=type(exc).__name__,
                diagnostics=getattr(exc, "diagnostics", None),
            )
            (out / "ERROR.txt").write_text(traceback.format_exc())
            tb, qp = exc.__traceback__, None
            while tb:
                frame = tb.tb_frame
                if frame.f_code.co_name == "solve" and frame.f_globals.get("__name__", "").endswith(
                    "osqp_motion"
                ):
                    qp = frame.f_locals
                tb = tb.tb_next
            if qp is not None:
                a, g = qp["a"], qp["g"]
                sparse.save_npz(out / "A.npz", a)
                sparse.save_npz(out / "G.npz", g)
                np.savez_compressed(
                    out / "QP_INPUT.npz",
                    b=qp["b"],
                    lower=qp["lower"],
                    primal_candidate=qp["u"],
                    multipliers_candidate=qp["lam"],
                )
                dense = g.toarray()  # small diagnostic SVD only; never motion production
                norm = np.linalg.norm(dense, axis=1)
                cosine = (dense / norm[:, None]) @ (dense / norm[:, None]).T
                np.fill_diagonal(cosine, 0)
                record["small_component_audit"] = dict(
                    dimensions=list(g.shape),
                    A_nnz=a.nnz,
                    G_nnz=g.nnz,
                    numerical_rank=int(np.linalg.matrix_rank(dense)),
                    rows_with_near_parallel_partner=int(np.sum(np.max(cosine, axis=1) > 1 - 1e-8)),
                    max_parallel_cosine=float(cosine.max()),
                    initial_empty_accepted_cache=not bool(tracker.translation.cache.particle_ids),
                )
        write(out / "DIAGNOSTIC.json", record)
        print(json.dumps(dict(case=summary["case"], status=record["status"], output=str(out))), flush=True)


if __name__ == "__main__":
    main()

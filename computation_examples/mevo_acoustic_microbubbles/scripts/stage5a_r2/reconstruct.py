"""Re-execute old rejected intervals read-only, retaining row/anchor metadata."""

from pathlib import Path
import json
import numpy as np
from scipy import sparse
from scripts.stage4l.dense_probe import load_case, sha, write
from scripts.stage5a_r2.seal import ROOT, BASE
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.bubbles.constrained_motion import PredictorCorrector
from mevo_mb.bubbles.translation_backends import _gap_vectors


def main():
    out = BASE / "diagnosis/reconstruction_v1"
    out.mkdir(exist_ok=False)
    reports = []
    old = ROOT / "results/stage5/motion_solver"
    for fixture in sorted((BASE / "fixtures").iterdir()):
        if not fixture.is_dir():
            continue
        prior = json.loads((old / "diagnostics" / fixture.name / "DIAGNOSTIC.json").read_text())
        original = Path(prior["source_run"])
        case = prior.get("case", "B2")
        cfg = json.loads((original / "CONFIG.json").read_text())
        _, field, flow, _, _ = load_case(Path("/workspace/mevo_acoustic_microbubbles"), case)
        tracker = ParticleTracker(cfg, field, flow)
        step = PredictorCorrector(tracker.translation, tracker.cylinder, tracker.evaluate)
        with np.load(original / "LAST_ACCEPTED_STATE.npz") as h:
            mask = h["active"]
            x = h["positions_m"][mask]
            r = h["radii_m"][mask]
            ids = h["particle_id"][mask]
            t = float(h["time_s"])
        if "failure_replay" in fixture.name:
            import h5py

            with h5py.File(original / "trajectories.h5") as h:
                xp = h["positions_m"][-2][mask]
                tp = float(h["time_s"][-2])
            _, _, cache, _ = step.trial(xp, r, ids, tp, t - tp)
            tracker.translation.cache = cache
            intervals = [prior["attempt_dt_s"]]
        else:
            rejections = [json.loads(s) for s in (original / "REJECTED.jsonl").read_text().splitlines()]
            intervals = [z["end_s"] - z["time_s"] for z in rejections] + [prior["dt_s"]]
        folder = out / fixture.name
        folder.mkdir()
        for index, dt in enumerate(intervals):
            row = dict(
                fixture=fixture.name, case=case, time_s=t, dt_s=dt, original_run=str(original), index=index
            )
            try:
                step.trial(x, r, ids, t, dt)
                row["status"] = "NOT_REPRODUCED"
            except Exception as exc:
                row.update(status=type(exc).__name__, diagnostics=getattr(exc, "diagnostics", {}))
                frames = {}
                tb = exc.__traceback__
                while tb:
                    frames[tb.tb_frame.f_code.co_name] = tb.tb_frame.f_locals.copy()
                    tb = tb.tb_next
                tr = frames.get("trial", {})
                co = frames.get("component", {})
                cg = frames.get("constrained", {})
                if "xp" in tr:
                    pairs = np.unique(np.vstack([tr["c0"].pairs, tr["cp"].pairs]), axis=0)
                    h0, n0 = _gap_vectors(tr["x"], pairs, r[0])
                    hp, npred = _gap_vectors(tr["xp"], pairs, r[0])
                    hc, nc = _gap_vectors(tr.get("xc", tr["xp"]), pairs, r[0])
                    denom = 0.002 * r[0] + 0.02 * np.maximum(h0, 0.001 * r[0])
                    pe = abs(hc - hp) / denom
                    winner = int(np.argmax(pe)) if len(pe) else None
                    row["responsible_pair"] = (
                        dict(
                            ids=ids[pairs[winner]].tolist(),
                            error_ratio=float(pe[winner]),
                            h_start_m=float(h0[winner]),
                            h_predict_m=float(hp[winner]),
                            h_correct_m=float(hc[winner]),
                            n_start=n0[winner].tolist(),
                            n_predict=npred[winner].tolist(),
                            n_correct=nc[winner].tolist(),
                            coefficient_start=float(
                                tracker.translation.coefficient(h0[winner : winner + 1], r[0])[0]
                            ),
                            coefficient_predict=float(
                                tracker.translation.coefficient(hp[winner : winner + 1], r[0])[0]
                            ),
                        )
                        if winner is not None
                        else None
                    )
                    row["switches"] = {
                        key: int(np.sum((h0 <= value) != (hp <= value)))
                        for key, value in [
                            ("contact_tolerance", 1e-11),
                            ("floor", 0.001 * r[0]),
                            ("activation", 0.1 * r[0]),
                        ]
                    }
                    row["position_responsible_id"] = int(
                        ids[np.argmax(np.linalg.norm(tr.get("xc", tr["xp"]) - tr["xp"], axis=1))]
                    )
                if index == len(intervals) - 1 and "aa" in co:
                    a, g = co["aa"], co["gg"]
                    saved_a = sparse.load_npz(fixture / "A.npz")
                    saved_g = sparse.load_npz(fixture / "G.npz")
                    row["reconstruction"] = dict(
                        A_exact=bool((a != saved_a).nnz == 0),
                        G_exact=bool((g != saved_g).nnz == 0),
                        A_max_difference=float(np.max(abs((a - saved_a).data), initial=0)),
                        G_max_difference=float(np.max(abs((g - saved_g).data), initial=0)),
                        velocity_scale_m_s=cg["vs"],
                        resistance_scale_N_s_m=cg["result"].scale,
                        force_scale_N=cg["vs"] * cg["result"].scale,
                        component_ids=ids[co["ix"]].tolist(),
                        full_row_indices=co["rr"].tolist(),
                        full_column_indices=co["cols"].tolist(),
                        warm_start="reconstructed preceding accepted trial"
                        if "failure_replay" in fixture.name
                        else "empty accepted cache",
                    )
                    lower = cg["lower"][co["rr"]] / cg["vs"]
                    up = tr["up"].ravel()[co["cols"]] / cg["vs"]
                    row["predictor_feasibility"] = dict(
                        min_slack=float(np.min(saved_g @ up - lower)),
                        zero_min_slack=float(np.min(-lower)),
                        positive_lower_rows=int(np.sum(lower > 0)),
                        start_min_slack=float(
                            min((saved_g @ up - lower)[[k[0] == "start" for k in co["keys"]]], default=np.inf)
                        ),
                        predictor_min_slack=float(
                            min(
                                (saved_g @ up - lower)[[k[0] == "predictor" for k in co["keys"]]],
                                default=np.inf,
                            )
                        ),
                    )
                    np.savez_compressed(
                        folder / "RECONSTRUCTED_CONTEXT.npz",
                        x=x,
                        xp=tr["xp"],
                        r=r,
                        ids=ids,
                        up=tr["up"],
                        predictor_normalized_component=up,
                        lower=lower,
                        cols=co["cols"],
                        rows=co["rr"],
                    )
                    keys = co["keys"]
                    meta = []
                    id_to_local = {int(k): i for i, k in enumerate(ids)}
                    for j, key in enumerate(keys):
                        stage, kind, *particle = key
                        anchor = tr["x"] if stage == "start" else tr["xp"]
                        locals_ = [id_to_local[k] for k in particle]
                        if kind == "pair":
                            v = anchor[locals_[1]] - anchor[locals_[0]]
                            norm = v / np.linalg.norm(v)
                            gap = np.linalg.norm(v) - 2 * r[0]
                            feature = None
                        else:
                            wall = tr["w0"] if stage == "start" else tr["wp"]
                            gap = wall[0][locals_[0]]
                            norm = wall[1][locals_[0]]
                            import pyvista as pv

                            surface = pv.read(cfg["geometry"]["mesh_surface"])
                            feature = int(surface.find_closest_cell(anchor[locals_[0]]))
                        meta.append(
                            dict(
                                original_row=j,
                                stage=stage,
                                canonical_contact=[kind, *particle],
                                normal=norm.tolist(),
                                anchor_positions_m=anchor[locals_].tolist(),
                                signed_gap_m=float(gap),
                                normalized_lower=float(lower[j]),
                                wall_exterior_cell_id=feature,
                                row_scaling=1.0,
                                working_set_row_mapping="native solver internal, not exposed",
                            )
                        )
                    write(folder / "ROW_METADATA.json", meta)
                    gd = saved_g.toarray()
                    s = np.linalg.svd(gd, compute_uv=False)
                    th = max(gd.shape) * np.finfo(float).eps * s[0]
                    cos = gd @ gd.T / np.outer(np.linalg.norm(gd, axis=1), np.linalg.norm(gd, axis=1))
                    np.fill_diagonal(cos, 0)
                    near = np.argwhere(np.triu(cos > 1 - 1e-8, 1))
                    row["spectrum"] = dict(
                        scaling="original saved normalized G; no row renormalization",
                        threshold=th,
                        singular_values=s.tolist(),
                        rank=int(np.sum(s > th)),
                        A_eigenvalues=np.linalg.eigvalsh(saved_a.toarray()).tolist(),
                        near_rows=[
                            dict(
                                i=int(i),
                                j=int(j),
                                cosine=float(cos[i, j]),
                                lower_difference=float(lower[i] - lower[j]),
                                row_difference_l2=float(np.linalg.norm(gd[i] - gd[j])),
                            )
                            for i, j in near
                        ],
                    )
            reports.append(row)
            write(out / "REPLAY.json", reports)
            print(fixture.name, index, row["status"], flush=True)
    write(
        out / "IDENTITY.json",
        dict(source_sha256=sha(Path(__file__)), fixture_manifest_sha256=sha(BASE / "fixtures/MANIFEST.json")),
    )


if __name__ == "__main__":
    main()

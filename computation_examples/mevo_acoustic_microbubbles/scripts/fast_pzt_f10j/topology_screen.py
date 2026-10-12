"""Historical data only; unchanged 2.1R graph and 2 ms persistence."""

import json
from collections import Counter
import h5py
import numpy as np
from mevo_mb.fast_sim.morphodynamics import topology
from .common import ROOT, BASE, write, js, sha


def main():
    rows = []
    skipped = []
    stages = ["f5", "f7", "f8"] + ["f10" + c for c in "abcdefghi"]
    for stage in stages:
        for path in sorted((ROOT / "results/fast_pzt" / stage).rglob("trajectory.h5")):
            rel = str(path.relative_to(ROOT))
            if any(t in rel for t in ["_serial_", "performance_v", "benchmark", "batch_v", "exclusive_v"]):
                skipped.append(
                    dict(
                        path=rel,
                        reason="replicated throughput or timing workload; not a new scientific topology scenario",
                    )
                )
                continue
            with h5py.File(path) as f:
                active_key = "active" if "active" in f else "trajectory/active"
                required = [
                    "positions_m",
                    "time_s",
                    active_key,
                    "particle_id",
                    "metrics/labels",
                    "edge_events/time_s",
                    "edge_events/batch_i_j_linked",
                ]
                missing = [name for name in required if name not in f]
                if missing:
                    skipped.append(
                        dict(
                            path=rel,
                            reason="Incomplete schema for fixed-ID accepted-ledger reconstruction",
                            missing=missing,
                        )
                    )
                    continue
                if "edge_events/time_s" not in f or f["positions_m"].shape[1] != 1:
                    skipped.append(dict(path=rel, reason="no complete accepted-step edge ledger or not B=1"))
                    continue
                times = f["time_s"][:]
                x = f["positions_m"][0, 0]
                ids = f["particle_id"][:]
                active = f[active_key][:]
                if not active.all():
                    skipped.append(
                        dict(path=rel, reason="active IDs change; not a complete fixed-ID topology window")
                    )
                    continue
                config = json.loads(f.attrs.get("config_json", "{}"))
                sp = path.parent / "SUMMARY.json"
                summary = js(sp) if sp.exists() else {}
                et = f["edge_events/time_s"][:]
                ei = f["edge_events/batch_i_j_linked"][:]
                try:
                    result, tt, labels = topology(
                        x, times, f["metrics/labels"][:, 0], et, ei, float(times[-1])
                    )
                    status = "ACCEPTED_LEDGER_RECONSTRUCTION_PASS"
                except (AssertionError, ValueError) as exc:
                    rows.append(dict(path=rel, status="SCHEMA_OR_LEDGER_NOT_RECONSTRUCTED", reason=str(exc)))
                    continue
                name = summary.get(
                    "numerical_mode", f.attrs.get("contact_model", config.get("mobility", "legacy M0"))
                )
                quality = summary.get("numerical_quality", f.attrs.get("numerical_quality", "UNSPECIFIED"))
                b2 = "B2_GEOM" in str(name)
                counts = np.array([len(np.unique(lab)) for lab in labels])
                rows.append(
                    dict(
                        path=rel,
                        sha256=sha(path),
                        stage=stage,
                        frames=len(times),
                        duration_s=float(times[-1]),
                        saved_cadence_max_s=float(np.max(np.diff(times))) if len(times) > 1 else 0.0,
                        n=len(ids),
                        particle_ids=ids,
                        model=name,
                        quality=str(quality),
                        status=status,
                        is_same_geometric_B2=b2,
                        complete=summary.get("status", str(f.attrs.get("complete", "unspecified"))),
                        component_count_range=[int(counts.min()), int(counts.max())],
                        topology=result,
                        compatibility="Same geometric B2, frozen-step resistance; no persistent events observed"
                        if b2
                        else "Historical geometry only; contact/mobility/controls differ. Not B2 numerical truth.",
                    )
                )
                if stage == "f5" and path.parent.name == "main":
                    first = next(
                        e
                        for e in result["events"]
                        if e["kind"] == "split"
                        and e["persistence_class"] == "persistent"
                        and e["time_s"] > 0.3
                    )
                    # A stored accepted frame strictly before the actual historical split.
                    index = max(0, int(np.searchsorted(times, first["time_s"] - 0.001)) - 1)
                    snapshot = dict(
                        positions_m=f["positions_m"][index, 0],
                        particle_ids=ids,
                        radii_m=f["radii_m"][:],
                        active=f["active"][index, 0],
                        wall_gap_m=f["wall_gap_m"][index, 0],
                        absolute_time_s=times[index],
                    )
                    np.savez_compressed(BASE / "analysis/HISTORICAL_F5_CANDIDATE.npz", **snapshot)
                    y = snapshot["positions_m"]
                    d = np.linalg.norm(y[:, None] - y[None, :], axis=-1)[np.triu_indices(len(y), 1)]
                    candidate = dict(
                        status="HISTORICAL_INITIAL_STATE_CANDIDATE_ONLY",
                        source=rel,
                        source_sha256=sha(path),
                        saved_frame=index,
                        absolute_start_time_s=float(times[index]),
                        historical_event=first,
                        proposed_window_s=[float(times[index]), float(first["time_s"] + 0.003)],
                        controls=config["events"],
                        applied_control="All PZT off; secondary attraction off; prescribed steady CFD on",
                        cfd="results/fast_pzt/f5/inputs/result_121.vtu",
                        wall="results/fast_pzt/f5/inputs/f4/prepared_geometry.npz",
                        input_manifest_sha256=sha(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
                        source_contact_model="F5 M0: local lubrication and positional correction; NOT B2",
                        min_wall_gap_m=float(np.min(snapshot["wall_gap_m"])),
                        max_overlap_m=float(max(0, 2 * 1.78e-6 - d.min())),
                        all_active=bool(snapshot["active"].all()),
                        finite_positions=bool(np.isfinite(y).all()),
                        persistence_s=0.002,
                        center_distance_threshold_over_R=2.1,
                        metrics=[
                            "fixed-ID component membership",
                            "persistent splits and merges with 2 ms dwell",
                            "right-censored events kept separate",
                            "position, shape, contact and wall budgets",
                        ],
                        prerequisites=[
                            "Separately preregister B2 restart semantics, coupling to near wall, quality and dt controls",
                            "Recheck full B2 first-step and wall budgets; F5 split time is not a reference target",
                        ],
                        B2_validation_executed=False,
                        training_ready=False,
                    )
                    write(BASE / "analysis/TOPOLOGY_CANDIDATE.json", candidate)
    summary = {
        s: dict(
            trajectories=sum(r.get("stage") == s for r in rows),
            persistent_splits=sum(
                r["topology"]["counts"]["split_persistent"] for r in rows if r.get("stage") == s
            ),
            persistent_merges=sum(
                r["topology"]["counts"]["merge_persistent"] for r in rows if r.get("stage") == s
            ),
        )
        for s in stages
    }
    write(
        BASE / "analysis/TOPOLOGY_SCREEN.json",
        dict(
            rows=rows,
            skipped=skipped,
            summary=summary,
            counts_are_correlated_historical_runs_not_independent_samples=True,
            status="PERSISTENT TOPOLOGY VALIDATION NOT RUN",
            new_trajectories=0,
            persistence_s=0.002,
            source_sha256=sha(__file__),
            reconstruction_status_counts=dict(Counter(r["status"] for r in rows)),
        ),
    )
    print(summary)


if __name__ == "__main__":
    main()

"""All-ID, actual-state analysis of the first B2 resolution-target crossing.

The six replays have two distinct real initial states. Within each lineage the
initial state is identical at all three levels; cross-lineage differences are
reported explicitly and never treated as same-initial-state solver errors.
"""

import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, read, write, sha
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r3.lineage import sample_all, metrics, summarize


def comparison(name, left, right, clock, group):
    a, b = [sample_all(p, clock) for p in (left, right)]
    assert np.array_equal(a["ids"], b["ids"])
    assert np.array_equal(a["radii"], b["radii"])
    groups, arrays = {}, {}
    for label, mask in [
        ("all_active", np.ones(len(a["ids"]), bool)),
        ("original_group", np.isin(a["ids"], group)),
        ("complement", ~np.isin(a["ids"], group)),
    ]:
        errors = metrics(a["positions"][:, mask], b["positions"][:, mask])
        groups[label] = dict(count=int(mask.sum()), errors_m={})
        for mode, values in errors.items():
            j, i = np.unravel_index(values.argmax(), values.shape)
            groups[label]["errors_m"][mode] = dict(
                initial=summarize(values[0]),
                endpoint=summarize(values[-1]),
                series=summarize(values),
                worst_id=int(a["ids"][mask][i]),
                worst_requested_time_s=float(clock[j]),
                worst_actual_times_s=[float(a["time"][j]), float(b["time"][j])],
            )
            arrays[f"{label}_{mode}_errors_m"] = values
    output = BASE / "comparisons" / f"B2_early_{name}.npz"
    np.savez_compressed(
        output,
        clock_s=clock,
        particle_id=a["ids"],
        radii_m=a["radii"],
        left_positions_m=a["positions"],
        right_positions_m=b["positions"],
        left_actual_time_s=a["time"],
        right_actual_time_s=b["time"],
        left_sampling_bound_m=a["uncertainty_m"],
        right_sampling_bound_m=b["uncertainty_m"],
        **arrays,
    )
    return dict(
        name=name,
        paths=[str(left), str(right)],
        sha256=[sha(left), sha(right)],
        all_active=len(a["ids"]),
        identical_local_initial_state=bool(np.array_equal(a["positions"][0], b["positions"][0])),
        groups=groups,
        endpoint_sampling_bound_m=float(a["uncertainty_m"][-1] + b["uncertainty_m"][-1]),
        max_sampling_bound_m=float(np.max(a["uncertainty_m"] + b["uncertainty_m"])),
        samples=str(output),
        no_coordinate_interpolation=True,
        accuracy_status="NOT_EVALUABLE",
    )


def contact_events(path, selected_ids=(82, 211)):
    """Actual accepted pair gaps and contact-band transitions, not force inference."""
    import pyvista as pv
    from scripts.stage5a_r2.motion_visuals import wall_path

    wall_source = wall_path()
    wall = pv.read(wall_source)
    with h5py.File(path) as h:
        times = h["time_s"][:]
        take = (times >= 3.00003 - 1e-14) & (times <= 3.00008 + 1e-14)
        x = h["positions_m"][:][take]
        times = times[take]
        mask = h["active"][0]
        ids, r = h["particle_id"][:][mask], h["radii_m"][:][mask]
        x = x[:, mask]
    events, wall_events, data = {}, {}, {"time_s": times, "particle_id": ids}
    for pid in selected_ids:
        i = int(np.flatnonzero(ids == pid)[0])
        gaps = np.linalg.norm(x - x[:, i, None], axis=2) - r[None, :] - r[i]
        gaps[:, i] = np.inf
        nearest = gaps.argmin(1)
        minimum = gaps[np.arange(len(times)), nearest]
        count = (gaps <= 1e-11).sum(1)
        changes = np.flatnonzero(np.r_[True, (nearest[1:] != nearest[:-1]) | (count[1:] != count[:-1])])
        events[str(pid)] = [
            dict(
                time_s=float(times[j]),
                nearest_particle_id=int(ids[nearest[j]]),
                raw_pair_gap_m=float(minimum[j]),
                pairs_within_positive_contact_tolerance=int(count[j]),
            )
            for j in changes
        ]
        data[f"id_{pid}_minimum_pair_gap_m"] = minimum
        data[f"id_{pid}_contact_band_count"] = count
        data[f"id_{pid}_nearest_id"] = ids[nearest]
        data[f"id_{pid}_position_m"] = x[:, i]
        cells, closest = wall.find_closest_cell(x[:, i], return_closest_point=True)
        displacement = x[:, i] - closest
        distance = np.linalg.norm(displacement, axis=1)
        normals = displacement / distance[:, None]
        jumps = np.linalg.norm(np.diff(normals, axis=0), axis=1)
        wall_events[str(pid)] = dict(
            minimum_unsigned_surface_gap_m=float(np.min(distance - r[i])),
            face_switch_count=int(np.sum(cells[1:] != cells[:-1])),
            maximum_nearest_direction_jump=float(jumps.max()) if len(jumps) else 0.0,
            face_ids=np.unique(cells).tolist(),
            semantics="Read-only finite WALL nearest feature at accepted coordinates; unsigned distance is not a signed/continuous safety certificate",
        )
        data[f"id_{pid}_wall_face_id"] = cells
        data[f"id_{pid}_unsigned_wall_gap_m"] = distance - r[i]
        data[f"id_{pid}_nearest_wall_direction"] = normals
    accepted = read(path.parent / "ACCEPTED.jsonl")
    rejected = read(path.parent / "REJECTED.jsonl")

    def interval(row):
        return 3.00003 - 1e-14 <= row["time_s"] <= 3.00008 + 1e-14

    relevant = [r for r in accepted if interval(r)]
    return dict(
        run=str(path.parent),
        accepted_steps=len(relevant),
        S0_steps=sum(r["substep_method"] == "S0" for r in relevant),
        fallback_calls=sum(r["fallbacks"] for r in relevant),
        rejection_records=[r for r in rejected if interval(r)],
        accepted_contact_band_transitions=events,
        wall_observations=wall_events,
        wall_source=str(wall_source),
        wall_source_sha256=sha(wall_source),
        caveat="Proximity-band changes are geometric observations, not an assertion that a QP multiplier became active; WALL feature/fallback causality requires independent records",
    ), data


def main():
    evidence = read(BASE / "diagnosis/EARLY_PRECISION.json")
    paths = {(r["kind"], r["level"]): local(r["run"]) / "trajectories.h5" for r in evidence["rows"]}
    assert all(r["summary"]["status"] == "COMPLETE" for r in evidence["rows"])
    group = np.load(R2 / "comparisons/B2_0.001_v4a_legacy_reference.npz")["particle_ids"]
    # These are common real output boundaries, not invented dense event frames.
    clock = np.array([3.00003, 3.00004, 3.00005])
    comparisons, orders = [], {}
    for kind in ("nominal", "refined"):
        for low, high in ((0, 1), (1, 2), (0, 2)):
            comparisons.append(
                comparison(f"{kind}_L{low}_L{high}", paths[kind, low], paths[kind, high], clock, group)
            )
        e0, e1 = [r["groups"]["all_active"]["errors_m"]["raw"]["endpoint"]["rms"] for r in comparisons[-3:-1]]
        orders[kind] = float(np.log2(e0 / e1)) if min(e0, e1) > 0 else None
    for level in range(3):
        comparisons.append(
            comparison(
                f"cross_initial_L{level}", paths["nominal", level], paths["refined", level], clock, group
            )
        )
    ledger = read(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")["ledger"]
    events = []
    for kind in ("nominal", "refined"):
        entry = next(
            r for r in ledger if r["case"] == "B2" and r["kind"] == kind and r["duration_s"] == 0.005
        )
        formal = local(entry["run"]) / "trajectories.h5"
        comparisons.append(comparison(f"{kind}_formal_vs_local_L2", formal, paths[kind, 2], clock, group))
        report, arrays = contact_events(formal)
        output = BASE / "comparisons" / f"B2_early_{kind}_contact_events.npz"
        np.savez_compressed(output, **arrays)
        report["arrays"] = str(output)
        events.append(report)
    write(
        BASE / "final/B2_EARLY_PRECISION_AUDIT.json",
        dict(
            status="NOT_EVALUABLE",
            source_sha256=sha(BASE / "diagnosis/EARLY_PRECISION.json"),
            runs=evidence["rows"],
            comparisons=comparisons,
            observed_endpoint_rms_orders=orders,
            formal_event_correlations=events,
            hypotheses={
                "S0_or_fallback_caused_first_crossing": {
                    "status": "REFUTED"
                    if all(e["S0_steps"] == e["fallback_calls"] == 0 for e in events)
                    else "UNRESOLVED",
                    "scope": "Actual +30..80 us accepted histories only",
                },
                "near_contact_rearrangement_correlates_with_divergence": {
                    "status": "CONFIRMED",
                    "scope": "Geometric contact-band and nearest-neighbor changes occur in both actual histories; correlation does not prove a unique cause",
                },
                "one_lineage_is_a_validated_global_truth": {"status": "UNRESOLVED"},
            },
            interpretation="Local three-level agreement and contact-band correlation cannot establish a global error bound or select one numerical lineage as truth",
            all_active_global_error_bound="NOT_ESTABLISHED",
            charged_experiment_wall_s=evidence["charged_experiment_wall_s"],
        ),
    )
    print({r["name"]: r["groups"]["all_active"]["errors_m"]["raw"]["endpoint"] for r in comparisons})


if __name__ == "__main__":
    main()

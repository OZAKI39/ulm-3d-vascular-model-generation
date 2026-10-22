"""Every accepted state, independent pair edges and unchanged persistent topology."""

from collections import defaultdict
from time import perf_counter
import h5py
import numpy as np
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.f10g_contact import shape, compare
from mevo_mb.fast_sim.morphodynamics import topology
from .common import ROOT, BASE, TRIALS, js, write, identity, read_steps, gz


def components(linked, pairs, n=32):
    a = np.eye(n, dtype=bool)
    i, j = pairs.T
    a[i, j] = a[j, i] = linked
    lab = connected_components(a, directed=False)[1]
    members = [np.flatnonzero(lab == k).tolist() for k in np.unique(lab)]
    return lab, members


def rebuild(times, gaps, pairs, cutoff):
    linked = gaps <= cutoff
    result = []
    for frame, p in np.argwhere(linked[1:] != linked[:-1]):
        k = int(frame + 1)
        before = float(gaps[k - 1, p])
        after = float(gaps[k, p])
        slope = (after - before) / (times[k] - times[k - 1])
        result.append(
            dict(
                pair=pairs[p].tolist(),
                on=int(linked[k, p]),
                frame=k,
                time_s=float(times[k]),
                bracket_s=times[k - 1 : k + 1],
                surface_gap_before_m=before,
                surface_gap_after_m=after,
                signed_threshold_before_m=before - cutoff,
                signed_threshold_after_m=after - cutoff,
                secant_normal_speed_m_s=float(slope),
                interpolated_crossing_estimate_s=float(times[k - 1] + (cutoff - before) / slope),
                interpolation_label="INTERPOLATED_CROSSING_ESTIMATE",
            )
        )
    return result, linked


def matched(a, b, tol):
    d = [defaultdict(list), defaultdict(list)]
    for table, rows in zip(d, [a, b]):
        for row in rows:
            table[tuple(row["pair"])].append(row)
    matches = []
    unmatched = []
    for pair in sorted(set(d[0]) | set(d[1])):
        aa, bb = d[0][pair], d[1][pair]
        if [e["on"] for e in aa] != [e["on"] for e in bb]:
            unmatched.append(dict(pair=pair, DT10=aa, DT5=bb))
            continue
        for occurrence, (x, y) in enumerate(zip(aa, bb)):
            delta = abs(x["time_s"] - y["time_s"])
            matches.append(
                dict(
                    pair_id=pair,
                    occurrence=occurrence,
                    on=x["on"],
                    DT10=x,
                    DT5=y,
                    time_difference_s=delta,
                    timing_gate="PASS" if delta <= tol + 1e-14 else "FAIL",
                    repeated_switches_on_pair=len(aa) > 1,
                    global_order_inversion_with=[],
                )
            )
    inversions = []
    for i, x in enumerate(matches):
        for j, y in enumerate(matches[:i]):
            da = x["DT10"]["time_s"] - y["DT10"]["time_s"]
            db = x["DT5"]["time_s"] - y["DT5"]["time_s"]
            if da * db < 0:
                strict = min(abs(da), abs(db)) > tol + 1e-14
                inversions.append(
                    dict(
                        events=[j, i],
                        separation_DT10_s=da,
                        separation_DT5_s=db,
                        exceeds_original_order_tolerance=strict,
                    )
                )
                x["global_order_inversion_with"].append(j)
                y["global_order_inversion_with"].append(i)
    return dict(
        events=matches,
        unmatched_pair_histories=unmatched,
        edge_id_sequence_pass=not unmatched,
        max_event_difference_s=max((e["time_difference_s"] for e in matches), default=None),
        timing_fail_count=sum(e["timing_gate"] == "FAIL" for e in matches),
        cross_pair_inversions=inversions,
        cross_pair_order_pass=not any(e["exceeds_original_order_tolerance"] for e in inversions),
    )


def band_runs(times, signed, band):
    inside = abs(signed) <= band
    begins = np.flatnonzero(inside & ~np.r_[False, inside[:-1]])
    ends = np.flatnonzero(inside & ~np.r_[inside[1:], False])
    rows = [
        dict(
            first_accepted_time_s=float(times[a]),
            last_accepted_time_s=float(times[b]),
            consecutive_inside_duration_s=float(times[b] - times[a]),
            accepted_points=int(b - a + 1),
            left_censored=bool(a == 0),
            right_censored=bool(b == len(times) - 1),
        )
        for a, b in zip(begins, ends)
    ]
    return dict(
        band_m=band,
        runs=rows,
        total_consecutive_inside_duration_s=sum(r["consecutive_inside_duration_s"] for r in rows),
        scope="Durations only between consecutive actually accepted points inside the band; no interpolated scientific frames",
    )


def main():
    start = perf_counter()
    p = js(ROOT / "configs/fast_pzt/f10i/protocol.json")
    r = p["radius_m"]
    cut = 0.1 * r
    pairs = np.column_stack(np.triu_indices(32, 1))
    pi = int(np.flatnonzero(np.all(pairs == [7, 24], axis=1))[0])
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    series = {}
    snapshots = {}
    trace_audit = []
    for mode, path in TRIALS.items():
        times, x, steps = read_steps(path)
        assert len(steps) == len(times) - 1 and all(s["accepted"] and s["error_bits"] == 0 for s in steps)
        distance = np.linalg.norm(x[:, pairs[:, 1]] - x[:, pairs[:, 0]], axis=-1)
        gap = distance - 2 * r
        events, linked = rebuild(times, gap, pairs, cut)
        with h5py.File(path / "trajectory.h5") as h:
            np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
            actual = [
                (round(float(t), 12), int(i), int(j), int(on))
                for t, (_, i, j, on) in zip(h["edge_events/time_s"][:], h["edge_events/batch_i_j_linked"][:])
            ]
            expected = [(round(e["time_s"], 12), *e["pair"], e["on"]) for e in events]
            assert sorted(actual) == sorted(expected) and len(events) == 29
            for t, v in zip(h["time_s"][:], h["positions_m"][:, 0]):
                k = int(np.argmin(abs(times - t)))
                assert abs(times[k] - t) < 1e-13
                np.testing.assert_array_equal(x[k], v)
        del steps
        soft = gz(path / "SOFT_STEPS.json.gz")
        assert len(soft) == len(times) - 1
        velocity = np.array([z["result"]["corrected_velocity_m_s"] for z in soft])
        contact = np.array([z["result"]["contact_force_N"] for z in soft])
        reactions = np.array([z["result"]["pair_force_N"] for z in soft])
        residual = np.array([z["result"]["mobility_equation_residual_m"] for z in soft])
        selected = {int(np.argmin(abs(times - t))) for t in p["cpu_saved_endpoint_times_s"]}
        selected.add(int(np.argmax(residual)) + 1)
        for k in selected:
            assert k > 0
            s = soft[k - 1]
            snapshots[f"{mode}_{k:04d}"] = dict(
                mode=mode,
                accepted_endpoint_s=times[k],
                step=k,
                initial_positions_m=s["initial_positions_m"],
                archived_result=s["result"],
            )
        for e in events:
            k = e["frame"]
            idx = int(np.flatnonzero(np.all(pairs == e["pair"], axis=1))[0])
            i, j = e["pair"]
            normal = (x[k, j] - x[k, i]) / distance[k, idx]
            e.update(
                endpoint_normal_speed_m_s=float(np.dot(normal, velocity[k - 1, j] - velocity[k - 1, i])),
                pair_geometric_spring_force_N=float(reactions[k - 1, idx]),
                minimum_absolute_threshold_distance_m=float(np.min(abs(gap[:, idx] - cut))),
            )
        normal = (x[1:, 24] - x[1:, 7]) / distance[1:, pi, None]
        pair_diag = dict(
            mode=mode,
            times_s=times,
            center_distance_m=distance[:, pi],
            surface_gap_m=gap[:, pi],
            signed_edge_distance_m=gap[:, pi] - cut,
            accepted_interval_normal_speed_m_s=np.diff(gap[:, pi]) / np.diff(times),
            endpoint_velocity_times_s=times[1:],
            endpoint_normal_speed_m_s=np.sum(normal * (velocity[:, 24] - velocity[:, 7]), axis=1),
            pair_geometric_spring_force_N=reactions[:, pi],
            total_contact_force_on_ID7_N=contact[:, 7],
            total_contact_force_on_ID24_N=contact[:, 24],
            events=[e for e in events if e["pair"] == [7, 24]],
            bands=[band_runs(times, gap[:, pi] - cut, b) for b in p["near_threshold_bands_m"]],
            max_archived_solver_residual_m=float(residual.max()),
        )
        labs = []
        members = []
        shapes = []
        for at, v in zip(linked, x):
            lab, member = components(at, pairs)
            labs.append(lab)
            members.append(member)
            shapes.append(shape(v, r))
        et = np.array([e["time_s"] for e in events])
        ei = np.array([[0, *e["pair"], e["on"]] for e in events])
        topo = topology(x[0], times, np.array(labs), et, ei, times[-1], persistence_s=p["persistence_s"])[0]
        series[mode] = dict(
            times=times,
            x=x,
            gap=gap,
            linked=linked,
            events=events,
            pair_diagnostic=pair_diag,
            labels=labs,
            members=members,
            shapes=shapes,
            topology=topo,
            residual=residual,
        )
        np.savez_compressed(
            out / f"ACCEPTED_GEOMETRY_{mode}.npz",
            time_s=times,
            positions_m=x,
            particle_id=np.arange(32),
            pairs=pairs,
            center_distance_m=distance,
            surface_gap_m=gap,
            signed_edge_distance_m=gap - cut,
            linked=linked,
            labels=labs,
            solver_residual_m=residual,
        )
        trace_audit.append(
            dict(mode=mode, exact_original_events=True, accepted_frames=len(times), raw_events=len(events))
        )
        del soft
    comparison = matched(series["DT10"]["events"], series["DT5"]["events"], p["edge_time_tolerance_s"])
    assert comparison["edge_id_sequence_pass"] and len(comparison["events"]) == 29
    assert abs(comparison["max_event_difference_s"] - 130e-6) < 1e-14
    write(
        out / "EVENT_MATCHES.json",
        dict(
            **comparison,
            reconstruction=trace_audit,
            raw_edge_timing="RAW_EDGE_TIMING_FAIL",
            edge_identity="EDGE_ID_SEQUENCE_PASS",
            source=identity(),
        ),
    )
    write(
        out / "PAIR_7_24.json",
        dict(
            pair=[7, 24],
            contact_distance_m=2 * r,
            edge_surface_gap_m=cut,
            series=[s["pair_diagnostic"] for s in series.values()],
        ),
    )
    import gzip
    import json
    from .common import clean

    with gzip.open(out / "CPU_SAVED_STATES.json.gz", "wt") as f:
        json.dump(clean(snapshots), f, separators=(",", ":"), allow_nan=False)
    values = []
    a, b = series.values()
    for i, t in enumerate(a["times"]):
        j = int(np.argmin(abs(b["times"] - t)))
        assert abs(t - b["times"][j]) < 1e-13
        u, v = a["shapes"][i], b["shapes"][j]
        values.append(
            dict(
                time_s=t,
                **compare(a["x"][i], b["x"][j]),
                Rg_relative=abs(v["Rg_m"] - u["Rg_m"]) / u["Rg_m"],
                aspect_relative=abs(v["aspect_ratio"] - u["aspect_ratio"]) / u["aspect_ratio"],
                tensor_relative=np.linalg.norm(np.array(v["shape_tensor_m2"]) - u["shape_tensor_m2"])
                / np.linalg.norm(u["shape_tensor_m2"]),
                members_equal=a["members"][i] == b["members"][j],
            )
        )
    maxima = {
        k: max(v[k] for v in values)
        for k in [
            "centroid_norm_m",
            "max_fixed_ID_m",
            "Kabsch_RMSD_m",
            "Rg_relative",
            "aspect_relative",
            "tensor_relative",
        ]
    }
    assert all(v["members_equal"] for v in values)
    for mode, s in series.items():
        assert all(m == [list(range(32))] for m in s["members"])
        assert len(s["topology"]["events"]) == 0
    write(
        out / "GROUP_TOPOLOGY.json",
        dict(
            comparison=values,
            maxima=maxima,
            modes={
                k: dict(
                    time_s=s["times"],
                    members=s["members"],
                    component_count=[len(m) for m in s["members"]],
                    largest_fraction=[max(map(len, m)) / 32 for m in s["members"]],
                    Rg_m=[z["Rg_m"] for z in s["shapes"]],
                    aspect=[z["aspect_ratio"] for z in s["shapes"]],
                    topology=s["topology"],
                )
                for k, s in series.items()
            },
            status="GROUP_LEVEL_NUMERICAL_EVIDENCE",
            persistent_status="NO_PERSISTENT_SPLIT_OR_MERGE_OBSERVED",
            predictive_persistent_event_accuracy="PERSISTENT_TOPOLOGY_NOT_EVALUATED",
            right_censored_events=0,
            scope="All 1502 accepted frames include the initial state twice. No observed multi-component or merge/split event; this is not validation of unseen splitting cases.",
            training_ready=False,
        ),
    )
    thresholds = []
    for factor in p["edge_distance_factors"]:
        threshold = (factor - 2) * r if factor != 2.1 else cut
        outputs = []
        for mode, s in series.items():
            ev, linked = rebuild(s["times"], s["gap"], pairs, threshold)
            members = [components(v, pairs)[1] for v in linked]
            outputs.append(
                dict(
                    mode=mode,
                    events=ev,
                    event_count=len(ev),
                    unique_event_pairs=len({tuple(e["pair"]) for e in ev}),
                    all_single_component=all(len(m) == 1 for m in members),
                    minimum_largest_fraction=min(max(map(len, m)) / 32 for m in members),
                    changed_edge_state_entries_vs_original=int(np.count_nonzero(linked != s["linked"])),
                )
            )
        thresholds.append(
            dict(
                factor=factor,
                surface_gap_m=threshold,
                modes=outputs,
                comparison=matched(outputs[0]["events"], outputs[1]["events"], p["edge_time_tolerance_s"]),
            )
        )
    write(
        out / "THRESHOLD_SENSITIVITY.json",
        dict(
            primary_factor=2.1,
            thresholds=thresholds,
            scope="Preregistered offline graph-definition perturbation only; no changed trajectory, forces, physics or original admission",
        ),
    )
    write(
        out / "OFFLINE_ANALYSIS.json",
        dict(
            source=identity(),
            body_s=perf_counter() - start,
            saved_cpu_states=len(snapshots),
            new_physical_steps=0,
        ),
    )
    print(
        "29/29 independently reconstructed; max lag",
        comparison["max_event_difference_s"],
        "timing failures",
        comparison["timing_fail_count"],
    )
    print("Group maxima", maxima)
    print(
        "Thresholds",
        [
            (x["factor"], [y["event_count"] for y in x["modes"]], x["comparison"]["max_event_difference_s"])
            for x in thresholds
        ],
    )


if __name__ == "__main__":
    main()

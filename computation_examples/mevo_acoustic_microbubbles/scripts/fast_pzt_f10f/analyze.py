"""Read-only actual substep, common-prefix morphology and all-attempt cost analysis."""

import gzip
import json
from collections import defaultdict
import h5py
import numpy as np
from mevo_mb.fast_sim.morphodynamics import rigid_fit, topology
from .common import BASE, HISTORY, ROOT, js, sha, write, identity


def gz(path):
    with gzip.open(path, "rt") as f:
        return json.load(f)


def overlap(x, r=1.78e-6):
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return max(0, float(2 * r - d.min()))


def local_rounds(proposal, rounds=12, r=1.78e-6):
    x = proposal.copy()
    positions = []
    res = []
    corrections = []
    for _ in range(rounds):
        dx = np.zeros_like(x)
        for i in range(len(x)):
            for j in range(len(x)):
                if i == j:
                    continue
                v = x[i] - x[j]
                d = np.linalg.norm(v)
                if d < 2 * r:
                    dx[i] += 0.25 * (2 * r - d) * v / d
        x = x + dx
        positions.append(x.copy())
        res.append(overlap(x))
        corrections.append(float(np.linalg.norm(dx, axis=1).max()))
    return positions, res, corrections


def read(path):
    with h5py.File(path / "trajectory.h5") as h:
        out = {
            k: h[k][:]
            for k in ["time_s", "positions_m", "edge_events/time_s", "edge_events/batch_i_j_linked"]
        }
        for group in ["metrics", "morphology", "trajectory"]:
            for k in h[group]:
                out[group + "/" + k] = h[group + "/" + k][:]
    return out


def compare(path_a, path_b):
    a, b = read(path_a), read(path_b)
    nt = min(len(a["time_s"]), len(b["time_s"]))
    assert np.array_equal(a["time_s"][:nt], b["time_s"][:nt])
    times = a["time_s"][:nt]
    end = float(times[-1])
    base = dict(
        reference=str(path_a.relative_to(ROOT)),
        candidate=str(path_b.relative_to(ROOT)),
        accepted_common_end_s=end,
        saved_frames=nt,
        source_sha256=[sha(path_a / "trajectory.h5"), sha(path_b / "trajectory.h5")],
        paired_acoustic_delta="NOT_EVALUATED_NO_SAME_MODE_OFF_CONTROL",
    )
    if nt < 2:
        return base | dict(
            status="NOT_EVALUATED_NO_NONZERO_COMMON_ACCEPTED_PREFIX",
            reason="Initial coordinates alone do not validate evolved morphology",
        )
    x, y = a["positions_m"][:nt, 0], b["positions_m"][:nt, 0]
    center = np.linalg.norm(x.mean(1) - y.mean(1), axis=1)
    axes = np.max(np.abs(x.mean(1) - y.mean(1)), axis=0)
    idmax = np.linalg.norm(y - x, axis=-1).max(1)
    kab = np.array([rigid_fit(xx, yy)[1] for xx, yy in zip(x, y)])

    def ar(key):
        return a[key][:nt, 0], b[key][:nt, 0]

    rg1, rg2 = ar("metrics/global_rg_m")
    asp1, asp2 = ar("morphology/fixed_aspect_ratio")
    ten1, ten2 = ar("morphology/fixed_tensor_m2")
    tensor = np.linalg.norm(ten2 - ten1, axis=(1, 2)) / np.maximum(np.linalg.norm(ten1, axis=(1, 2)), 1e-300)
    large1, large2 = ar("metrics/largest_size")
    gap1, gap2 = ar("metrics/minimum_gap_m")
    diag1, diag2 = ar("trajectory/diag")
    pathdelta = np.abs(diag1[:, :, 4] - diag2[:, :, 4])
    pathlimit = 1e-9 + 0.02 * np.abs(diag1[:, :, 4])
    edges = []
    for data in [a, b]:
        mask = data["edge_events/time_s"] <= end + 1e-14
        edges.append((data["edge_events/batch_i_j_linked"][mask], data["edge_events/time_s"][mask]))
    bypair = []
    for ids, tt in edges:
        d = defaultdict(list)
        for (_, i, j, on), t in zip(ids, tt):
            d[(int(i), int(j))].append((int(on), float(t)))
        bypair.append(d)
    same_id = set(bypair[0]) == set(bypair[1])
    max_event_dt = 0.0
    paired = []
    missing = 0
    for key in set(bypair[0]) | set(bypair[1]):
        aa, bb = bypair[0][key], bypair[1][key]
        if [v[0] for v in aa] != [v[0] for v in bb]:
            same_id = False
            missing += 1
            continue
        for va, vb in zip(aa, bb):
            max_event_dt = max(max_event_dt, abs(va[1] - vb[1]))
            paired.append((va[1], vb[1]))
    order_ok = same_id
    for i in range(len(paired)):
        for j in range(i):
            da = paired[i][0] - paired[j][0]
            db = paired[i][1] - paired[j][1]
            if da * db < 0 and min(abs(da), abs(db)) > 1e-5 + 1e-14:
                order_ok = False
    tops = []
    for data, (ids, tt) in zip([a, b], edges):
        top, _, _ = topology(data["positions_m"][0, 0], times, data["metrics/labels"][:nt, 0], tt, ids, end)
        tops.append(top)

    def persistent(t):
        return [
            (e["kind"], e["particle_ids"], e["groups"], e["time_s"])
            for e in t["events"]
            if e["persistence_class"] == "persistent"
        ]

    pa, pb = persistent(tops[0]), persistent(tops[1])
    persistent_ok = len(pa) == len(pb) and all(
        u[:3] == v[:3] and abs(u[3] - v[3]) <= 1e-5 + 1e-14 for u, v in zip(pa, pb)
    )
    membership = np.array_equal(a["metrics/labels"][:nt], b["metrics/labels"][:nt])
    values = dict(
        centroid_max_m=float(center.max()),
        centroid_axis_max_m=axes,
        fixed_ID_max_m=float(idmax.max()),
        Kabsch_max_m=float(kab.max()),
        Rg_max_abs_m=float(np.abs(rg2 - rg1).max()),
        Rg_max_relative=float(np.max(np.abs(rg2 - rg1) / rg1)),
        aspect_max_abs=float(np.abs(asp2 - asp1).max()),
        aspect_max_relative=float(np.max(np.abs(asp2 - asp1) / asp1)),
        tensor_max_relative=float(tensor.max()),
        largest_fraction_max_abs=float(np.abs(large2 - large1).max() / 32),
        minimum_gap_max_abs_m=float(np.abs(gap2 - gap1).max()),
        contact_path_max_abs_m=float(pathdelta.max()),
        contact_path_max_gate_ratio=float((pathdelta / pathlimit).max()),
        event_count=[len(e[1]) for e in edges],
        edge_identity_per_pair_order=same_id,
        unmatched_pair_histories=missing,
        matched_edge_max_time_difference_s=max_event_dt,
        cross_pair_order=order_ok,
        persistent_events_equal=persistent_ok,
        component_membership_equal=membership,
        topology=tops,
    )
    gates = dict(
        centroid=values["centroid_max_m"] <= 1e-8,
        fixed_ID=values["fixed_ID_max_m"] <= 1.78e-7,
        Kabsch=values["Kabsch_max_m"] <= 8.9e-8,
        Rg=values["Rg_max_relative"] <= 0.01,
        aspect=values["aspect_max_relative"] <= 0.02,
        tensor=values["tensor_max_relative"] <= 0.02,
        largest_fraction=values["largest_fraction_max_abs"] == 0,
        minimum_gap=values["minimum_gap_max_abs_m"] <= 1.78e-8,
        contact_path=values["contact_path_max_gate_ratio"] <= 1,
        component_membership=membership,
        edge_events=same_id and order_ok and max_event_dt <= 1e-5 + 1e-14,
        persistent_events=(
            "NOT_EVALUATED_RIGHT_CENSORED"
            if any(e["persistence_class"] == "right_censored" for t in tops for e in t["events"])
            else persistent_ok
        ),
    )
    return base | dict(
        status="MORPHOLOGY_NOT_CONVERGED"
        if not all(v is True for v in gates.values())
        else "COMMON_PREFIX_ONLY_PASS",
        values=values,
        gates=gates,
        series=dict(
            time_s=times,
            centroid_difference_m=center,
            fixed_ID_max_m=idmax,
            Kabsch_RMSD_m=kab,
            Rg_reference_m=rg1,
            Rg_candidate_m=rg2,
            aspect_reference=asp1,
            aspect_candidate=asp2,
        ),
    )


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    data = js(BASE / "short_v1/RESULTS.json")
    rows = []
    local = []
    for row in data["rows"]:
        path = BASE / "short_v1" / (row["case"] + "_" + row["mode"])
        trace = gz(path / "CONTACT_TRACE.json.gz")
        ledger = gz(path / "SUBSTEP_LEDGER.json.gz")["steps"]
        last = trace["traces"][-1]
        dt = row["setting"]["dt_s"]
        iters = row["setting"]["contact_iterations"]
        assert not ledger[-1]["accepted"] and all(s["accepted"] for s in ledger[:-1])
        nsteps = len(ledger)
        rows.append(
            dict(
                case=row["case"],
                mode=row["mode"],
                first_rejection_s=ledger[-1]["attempted_end_s"],
                last_legal_s=ledger[-1]["start_s"],
                attempted_steps=nsteps,
                accepted_steps=nsteps - 1,
                final_predicted_nm=ledger[-1]["predicted_overlap_m"] * 1e9,
                first_predicted_nm=ledger[0]["predicted_overlap_m"] * 1e9,
                final_residual_nm=ledger[-1]["final_overlap_m"] * 1e9,
                final_rounds_nm=np.array(ledger[-1]["round_residual_m"]) * 1e9,
                final_passes_nm=np.array(ledger[-1]["round_max_correction_m"]) * 1e9,
                maximum_prediction_nm=max(s["predicted_overlap_m"] for s in ledger) * 1e9,
                maximum_pass_nm=max(max(s["round_max_correction_m"]) for s in ledger) * 1e9,
                maximum_accepted_residual_nm=max([s["final_overlap_m"] for s in ledger[:-1]] or [0]) * 1e9,
                maximum_wall_correction_nm=max(s["wall_correction_max_m"] for s in ledger) * 1e9,
                flow_queries=2 * nsteps,
                acoustic_updates=nsteps,
                nearfield_updates=nsteps,
                mobility_passes=2 * nsteps,
                pbd_passes=iters * nsteps,
                wall_passes=nsteps * (1 + iters // 3),
                integration_ms_per_10us=row["timing"]["integration_quality_wall_s"]
                / (nsteps * dt / 1e-5)
                * 1e3,
                timing=row["timing"],
                ledger_sha256=sha(path / "SUBSTEP_LEDGER.json.gz"),
            )
        )
        if row["mode"] == "B0":
            pred = np.array(last["predicted_positions_m"])
            pos, res, corr = local_rounds(pred)
            err = max(
                float(np.max(np.abs(x - np.array(y["positions_m"])))) for x, y in zip(pos, last["rounds"])
            )
            assert err < 1e-15
            local.append(
                dict(
                    case=row["case"],
                    scope="Same saved B0 rejected prediction; read-only Jacobi, no time advancement",
                    initial_prediction_nm=overlap(pred) * 1e9,
                    residual_nm=np.array(res) * 1e9,
                    max_pass_nm=np.array(corr) * 1e9,
                    max_numpy_gpu_position_difference_m=err,
                )
            )
    write(
        out / "CONTACT.json",
        dict(
            rows=rows,
            same_proposal=local,
            caution="First-failure depths occur at different physical times. Only O first-substep predictions share the exact input state.",
            O_common_initial_prediction_reduction_percent=[
                100
                * (
                    1
                    - next(r for r in rows if r["case"] == "O_1000_1102" and r["mode"] == m)[
                        "first_predicted_nm"
                    ]
                    / next(r for r in rows if r["case"] == "O_1000_1102" and r["mode"] == "B0")[
                        "first_predicted_nm"
                    ]
                )
                for m in ["SS2", "SS4"]
            ],
        ),
    )
    morph = []
    for case in ["prepare_1104", "O_1000_1102"]:
        for a, b in [("B0", "SS2"), ("B0", "SS4"), ("SS2", "SS4")]:
            morph.append(compare(BASE / "short_v1" / (case + "_" + a), BASE / "short_v1" / (case + "_" + b)))
    for m in ["SS2", "SS4"]:
        morph.append(
            compare(HISTORY / "extension_v1/O_1000_1102_C1", BASE / "short_v1" / ("O_1000_1102_" + m))
        )
    # Historical C1 agreement is a separate reference, never new F10F success.
    morph.append(compare(HISTORY / "contact_v1/O_1000_1102_C0", HISTORY / "extension_v1/O_1000_1102_C1"))
    write(out / "MORPHOLOGY.json", dict(comparisons=morph, training_ready=False))
    decision = dict(
        full_extension_admitted=False,
        SS12_admitted=False,
        batch_admitted=False,
        reason="All four SS2/SS4 scientific attempts failed original residual budget before original failure windows; no candidate passed the first contact admission gate.",
        formal_matched_performance="NOT_ADMITTED",
        further_seed_tuning=False,
        training_ready=False,
    )
    write(BASE / "STOP_DECISION.json", decision)
    process = [dict(label=p.name, **js(p)) for p in BASE.glob("*_PROCESS.json")]
    full = sum(p["full_process_wall_s"] for p in process)
    historic = js(HISTORY / "analysis/EFFECTIVE_YIELD.json")
    prior = historic["all_process_wall_interval_s"]
    yielddata = dict(
        new_complete_10ms_numerical_trajectories=0,
        new_independent_seeds=0,
        reused_independent_seed_families=2,
        complete_same_mode_OFF_S_O_triplets=0,
        scientific_attempts=6,
        scientific_contact_failures=6,
        smallsteps_attempts=4,
        smallsteps_failures=4,
        attempted_physical_steps=sum(r["attempted_steps"] for r in rows),
        accepted_physical_steps=sum(r["accepted_steps"] for r in rows),
        numerical_modes_tested=3,
        diagnostic_body_s=sum(r["timing"]["resident_body_s"] for r in rows),
        static_prepare_s=data["static_prepare_s"],
        cold_all_process_s=full,
        processes=process,
        actual_complete_trajectories_per_reserved_gpu_hour=0,
        complete_triplets_per_reserved_gpu_hour=0,
        seconds_per_new_qualified_sample=None,
        seconds_per_new_qualified_sample_status="UNDEFINED_ZERO_YIELD",
        cumulative_F10D_E_F_existing_complete_control_trajectories=3,
        cumulative_F10D_E_F_process_cost_interval_s=[v + full for v in prior],
        cumulative_F10D_E_F_existing_trajectories_per_hour_interval=[
            3 * 3600 / (v + full) for v in prior[::-1]
        ],
        production_speedup="NOT_EVALUATED_NO_ADMITTED_SMALLSTEPS_CANDIDATE",
        subkernel_time="NOT_EVALUATED: no formal timing experiment admitted; diagnostic CUDA trajectory spans include CPU gaps",
        count_scope="Actual attempted physical steps from per-step ledger, including first rejection. Counts follow unchanged fresh-physics call schedule; per-particle work is 32 times call counts where appropriate.",
        time_scope="Exclusive card reservation including process, diagnostics and output; not pure GPU kernel time or steady production. Engineering test costs retained.",
        training_ready=False,
    )
    write(out / "COST_AND_YIELD.json", yielddata)
    write(
        out / "SOURCE_MANIFEST.json",
        dict(
            analysis_source=identity(),
            inputs={
                str(p.relative_to(ROOT)): sha(p)
                for path in [BASE / "short_v1", HISTORY / "extension_v1", HISTORY / "contact_v1"]
                for p in path.rglob("*")
                if p.is_file()
            },
        ),
    )
    print("contact", [(r["case"], r["mode"], r["first_rejection_s"]) for r in rows])
    for m in morph:
        print(m["candidate"].split("/")[-1], m["status"], m.get("values", {}).get("fixed_ID_max_m"))


if __name__ == "__main__":
    main()

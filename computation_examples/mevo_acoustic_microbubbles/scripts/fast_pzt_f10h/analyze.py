"""Offline per-step resistance audit and fixed-ID comparisons on actual physical times."""

import gzip
import json
from collections import defaultdict
from time import perf_counter
import numpy as np
import h5py
from scipy import linalg
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.f10g_contact import resistance, balance, compare, shape
from mevo_mb.fast_sim.morphodynamics import topology
from .common import ROOT, BASE, js, write, sha, identity


def gz(p):
    with gzip.open(p, "rt") as f:
        return json.load(f)


def read_steps(path):
    if (path / "ALL_STEPS.json.gz").exists():
        rows = gz(path / "ALL_STEPS.json.gz")
        good = [s for s in rows if s["accepted"]]
        x = np.array([rows[0]["initial_positions_m"]] + [s["positions_after_step_m"] for s in good])
        times = np.array([0] + [s["attempted_end_s"] for s in good])
    else:
        rows = gz(path / "SOFT_STEPS.json.gz")
        good = [s for s in rows if s["accepted"]]
        x = np.array([rows[0]["initial_positions_m"]] + [s["actual_accepted_positions_m"] for s in good])
        times = np.r_[0, np.arange(1, len(good) + 1) * rows[0]["result"]["dt_s"]]
    return times, x, rows


def geometric_series(path):
    times, x, rows = read_steps(path)
    shapes = [shape(v, 1.78e-6) for v in x]
    labels = np.array(
        [
            connected_components(
                np.linalg.norm(v[:, None] - v[None, :], axis=-1) <= 2.1 * 1.78e-6, directed=False
            )[1]
            for v in x
        ]
    )
    with h5py.File(path / "trajectory.h5") as h:
        et = h["edge_events/time_s"][:]
        ids = h["edge_events/batch_i_j_linked"][:]
        mask = et <= times[-1] + 1e-14
        et = et[mask]
        ids = ids[mask]
        for t, v in zip(h["time_s"][:], h["positions_m"][:, 0]):
            at = int(np.argmin(abs(times - t)))
            assert abs(times[at] - t) < 1e-13
            np.testing.assert_array_equal(v, x[at])
    top, _, _ = topology(x[0], times, labels, et, ids, times[-1])
    return dict(times=times, x=x, shapes=shapes, labels=labels, edge_times=et, edge_ids=ids, topology=top)


def compare_paths(a, b):
    sa, sb = geometric_series(a), geometric_series(b)
    ta, tb = sa["times"], sb["times"]
    end = min(ta[-1], tb[-1])
    pairs = []
    for i, t in enumerate(ta):
        if t > end + 1e-14:
            break
        j = int(np.argmin(abs(tb - t)))
        if abs(tb[j] - t) < 1e-13:
            pairs.append((i, j))
    values = []
    for i, j in pairs:
        x, y = sa["x"][i], sb["x"][j]
        u, v = sa["shapes"][i], sb["shapes"][j]

        def largest(s):
            return sorted([m for m in s["members"] if len(m) == max(map(len, s["members"]))])

        values.append(
            dict(
                time_s=ta[i],
                **compare(x, y),
                Rg_relative=abs(v["Rg_m"] - u["Rg_m"]) / u["Rg_m"],
                aspect_relative=abs(v["aspect_ratio"] - u["aspect_ratio"]) / u["aspect_ratio"],
                tensor_relative=np.linalg.norm(np.array(v["shape_tensor_m2"]) - u["shape_tensor_m2"])
                / np.linalg.norm(u["shape_tensor_m2"]),
                largest_members_equal=largest(u) == largest(v),
                all_members_equal=u["members"] == v["members"],
                Rg_a_m=u["Rg_m"],
                Rg_b_m=v["Rg_m"],
                aspect_a=u["aspect_ratio"],
                aspect_b=v["aspect_ratio"],
            )
        )
    eventsets = []
    tops = []
    for s in [sa, sb]:
        mask = s["edge_times"] <= end + 1e-14
        et = s["edge_times"][mask]
        ids = s["edge_ids"][mask]
        d = defaultdict(list)
        for (_, i, j, on), t in zip(ids, et):
            d[(int(i), int(j))].append((int(on), float(t)))
        eventsets.append(d)
        maskt = s["times"] <= end + 1e-14
        tops.append(topology(s["x"][0], s["times"][maskt], s["labels"][maskt], et, ids, end)[0])
    unmatched = []
    maxdt = 0.0
    matched = []
    for pair in sorted(set(eventsets[0]) | set(eventsets[1])):
        u, v = [d[pair] for d in eventsets]
        if [z[0] for z in u] != [z[0] for z in v]:
            unmatched.append(dict(pair=pair, a=u, b=v))
        else:
            for aa, bb in zip(u, v):
                maxdt = max(maxdt, abs(aa[1] - bb[1]))
                matched.append([aa[1], bb[1]])
    order = True
    for i, (ta1, tb1) in enumerate(matched):
        for ta2, tb2 in matched[:i]:
            if (ta1 - ta2) * (tb1 - tb2) < 0 and min(abs(ta1 - ta2), abs(tb1 - tb2)) > 1e-5 + 1e-14:
                order = False
    persist = [[e for e in t["events"] if e["persistence_class"] == "persistent"] for t in tops]
    pequal = len(persist[0]) == len(persist[1]) and all(
        all(u[k] == v[k] for k in ["kind", "particle_ids", "groups"])
        and abs(u["time_s"] - v["time_s"]) <= 1e-5 + 1e-14
        for u, v in zip(*persist)
    )
    censored = any(e["persistence_class"] == "right_censored" for t in tops for e in t["events"])
    maximum = {
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
    threshold = js(ROOT / "configs/fast_pzt/f10h/protocol.json")["morphology"]
    gates = {
        k: bool(maximum[key] <= threshold[k])
        for k, key in [
            ("centroid_m", "centroid_norm_m"),
            ("fixed_ID_m", "max_fixed_ID_m"),
            ("Kabsch_m", "Kabsch_RMSD_m"),
            ("Rg_relative", "Rg_relative"),
            ("aspect_relative", "aspect_relative"),
            ("tensor_relative", "tensor_relative"),
        ]
    }
    gates.update(
        largest_members=all(v["largest_members_equal"] for v in values),
        all_members=all(v["all_members_equal"] for v in values),
        edge_events=not unmatched and order and maxdt <= 1e-5 + 1e-14,
        persistent_events="NOT_EVALUATED_RIGHT_CENSORED" if censored else pequal,
    )
    return dict(
        reference=str(a.relative_to(ROOT)),
        candidate=str(b.relative_to(ROOT)),
        source_sha256=[sha(a / "trajectory.h5"), sha(b / "trajectory.h5")],
        common_end_s=end,
        common_times=len(values),
        maxima=maximum,
        gates=gates,
        status="COMMON_PREFIX_PASS" if all(v is True for v in gates.values()) else "MORPHOLOGY_NOT_CONVERGED",
        series=values,
        event_counts=[sum(map(len, d.values())) for d in eventsets],
        unmatched_pair_histories=unmatched,
        matched_edge_max_time_difference_s=maxdt,
        cross_pair_order=order,
        topology=tops,
        scope="Same mode 10/5us only when explicitly selected; comparisons between laws are descriptive, not self-convergence. True physical times, no interpolation.",
    )


def audit_lub(path):
    start = perf_counter()
    times, x, steps = read_steps(path)
    r = 1.78e-6
    records = []
    for s in steps:
        xx = np.array(s["initial_positions_m"])
        free = np.array(s["free_velocity_m_s"])
        raw = np.array(s["local_velocity_m_s"])
        op = resistance(xx, r)
        ub = balance(op, free, raw)
        ref = linalg.solve(op["A"], free.ravel(), assume_a="pos").reshape(-1, 3)
        pred = np.array(s["predicted_positions_m"])
        assert np.max(abs(xx + s["dt_s"] * ub - pred)) < 1e-15
        residual = (op["A"] @ ub.ravel() - free.ravel()).reshape(-1, 3)
        pre = np.linalg.solve(op["blocks"], residual[..., None])[..., 0]
        scale = np.linalg.norm(free, axis=1).max()
        err = ub - ref
        i, j = op["pairs"].T
        normal = np.sum((ub[j] - ub[i]) * op["normals"], axis=1)
        nerr = np.sum((err[j] - err[i]) * op["normals"], axis=1)
        records.append(
            dict(
                time_s=s["attempted_end_s"],
                accepted=s["accepted"],
                normalized_preconditioned=float(np.linalg.norm(pre, axis=1).max() / scale),
                normalized_velocity_error=float(np.linalg.norm(err, axis=1).max() / scale),
                residual_max_m_s=float(np.linalg.norm(residual, axis=1).max()),
                residual_max_N=float(6 * np.pi * 0.001 * r * np.linalg.norm(residual, axis=1).max()),
                max_approach_speed_m_s=max(0, -float(normal.min(initial=0))),
                max_normal_error_m_s=float(np.max(abs(nerr), initial=0)),
                worst_normal_pair=op["pairs"][np.argmax(abs(nerr))] if len(nerr) else [],
                mean_velocity_m_s=ub.mean(0),
                floor_pair_count=int(op["floored"].sum()),
                accepted_displacement_max_m=float(
                    np.linalg.norm(np.array(s["positions_after_step_m"]) - xx, axis=1).max()
                ),
                velocity_error_per_ID_m_s=err,
                preconditioned_per_ID_m_s=pre,
            )
        )
    geom = geometric_series(path)
    ledger = gz(path / "SUBSTEP_LEDGER.json.gz")["steps"]
    good = [s for s in ledger if s["accepted"]]
    return dict(
        path=str(path.relative_to(ROOT)),
        accepted_end_s=times[-1],
        accepted_steps=len(good),
        attempted_steps=len(steps),
        max_accepted_overlap_m=max((s["final_overlap_m"] for s in good), default=0),
        max_single_pass_m=max((max(s["round_max_correction_m"]) for s in ledger), default=0),
        max_normalized_preconditioned=max(v["normalized_preconditioned"] for v in records),
        max_normalized_velocity_error=max(v["normalized_velocity_error"] for v in records),
        records=records,
        per_step_shapes=geom["shapes"],
        topology=geom["topology"],
        offline_analysis_s=perf_counter() - start,
    )


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    trials = []
    for label in ["short_v1", "long_v1", "dt5_v1"]:
        if (BASE / label / "RESULTS.json").exists():
            trials += js(BASE / label / "RESULTS.json")["rows"]
    audits = []
    for row in trials:
        audits.append(audit_lub(ROOT / row["path"]))
    write(out / "LUBRICATION_TRAJECTORIES.json", dict(rows=audits, trials=trials, source=identity()))
    comps = []
    for row in trials:
        if row["mode"].endswith("DT5"):
            ref = next(
                s
                for s in trials
                if s["case"] == row["case"]
                and s["mode"].endswith("DT10")
                and s["target_s"] == row["target_s"]
            )
            comps.append(
                compare_paths(ROOT / ref["path"], ROOT / row["path"])
                | dict(case=row["case"], comparison="SAME_LUB4_MODEL_DT10_DT5")
            )
    write(
        out / "MORPHOLOGY.json",
        dict(
            comparisons=comps,
            thresholds=js(ROOT / "configs/fast_pzt/f10h/protocol.json")["morphology"],
            training_ready=False,
        ),
    )
    b2_rows = []
    b2_trials = []
    for label in ["b2_short_v1", "b2_dt5_v1"]:
        if not (BASE / label / "RESULTS.json").exists():
            continue
        b2_trials += js(BASE / label / "RESULTS.json")["rows"]
        for trial in js(BASE / label / "RESULTS.json")["rows"]:
            path = ROOT / trial["path"]
            steps = gz(path / "SOFT_STEPS.json.gz")
            ledger = gz(path / "SUBSTEP_LEDGER.json.gz")["steps"]
            good = [s for s in steps if s["accepted"]]
            geom = geometric_series(path)
            false = []
            law = []
            torque = []
            iters = []
            for s in steps:
                v = s["result"]
                f = np.array(v["pair_force_N"])
                delta = np.array(v["geometric_overlap_m"])
                x = np.array(v["positions_m"])
                force = np.array(v["contact_force_N"])
                for pair, fr in zip(np.array(v["pairs"])[(delta == 0) & (f > 0)], f[(delta == 0) & (f > 0)]):
                    false.append(dict(step=s["step"], pair=pair, force_N=fr))
                law.append(float(np.max(abs(f - v["kc_N_m"] * delta))))
                torque.append(float(np.linalg.norm(np.cross(x - x.mean(0), force).sum(0))))
                iters.append(len(v["iterations"]))
            b2_rows.append(
                dict(
                    case=trial["case"],
                    mode=trial["mode"],
                    path=trial["path"],
                    accepted_steps=len(good),
                    accepted_end_s=geom["times"][-1],
                    max_accepted_overlap_m=max(max(s["result"]["geometric_overlap_m"]) for s in good),
                    max_single_iteration_m=max(s["result"]["max_single_iteration_m"] for s in steps),
                    max_equation_residual_m=max(s["result"]["mobility_equation_residual_m"] for s in steps),
                    max_spring_law_force_error_N=max(law),
                    false_reactions=false,
                    max_torque_N_m=max(torque),
                    max_net_force_N=max(
                        float(np.linalg.norm(s["result"]["net_contact_force_N"])) for s in steps
                    ),
                    max_contact_centroid_shift_m=max(
                        float(np.linalg.norm(s["result"]["centroid_contact_shift_m"])) for s in steps
                    ),
                    mean_outer_iterations=float(np.mean(iters)),
                    max_outer_iterations=max(iters),
                    step_wall_s=sum(s["full_step_s"] for s in steps),
                    field_and_read_s=sum(s["field_and_read_s"] for s in steps),
                    resistance_assembly_s=sum(s["resistance_assembly_s"] for s in steps),
                    contact_initialization_s=sum(s["result"]["initialization_s"] for s in steps),
                    contact_solver_s=sum(s["result"]["contact_s"] for s in steps),
                    commit_and_read_s=sum(s["commit_and_read_s"] for s in steps),
                    max_wall_correction_m=max(s["wall_correction_max_m"] for s in ledger),
                    topology=geom["topology"],
                    per_step_shapes=geom["shapes"],
                )
            )
    b2_comps = []
    if len(b2_trials) == 2:
        b2_comps.append(
            compare_paths(ROOT / b2_trials[0]["path"], ROOT / b2_trials[1]["path"])
            | dict(case="O_1000_1102", comparison="SAME_GEOMETRIC_B2_MODEL_DT10_DT5")
        )
    write(
        out / "B2_GEOMETRY_TRAJECTORIES.json",
        dict(rows=b2_rows, trials=b2_trials, comparisons=b2_comps, training_ready=False),
    )
    # Existing models are compared separately; these are not self-convergence errors.
    from scripts.fast_pzt_f10f.analyze import compare as historical_compare

    between = []
    if b2_trials:
        b2 = ROOT / b2_trials[0]["path"]
        lub = next(
            r
            for r in trials
            if r["case"] == "O_1000_1102" and r["mode"].endswith("DT10") and r["target_s"] == 0.005
        )
        for ref in [ROOT / lub["path"], ROOT / "results/fast_pzt/f10g/integration_v1/O_1000_1102_B2_K2"]:
            c = historical_compare(ref, b2)
            c["interpretation"] = (
                "Different contact/mobility algorithms; descriptive difference, not temporal convergence or physical truth"
            )
            between.append(c)
    write(out / "BETWEEN_MODELS.json", dict(comparisons=between))
    cost_rows = (
        js(BASE / "performance_v1/RESULTS.json")["rows"]
        if (BASE / "performance_v1/RESULTS.json").exists()
        else []
    )
    medians = {}
    for mode in sorted({s["mode"] for s in cost_rows}):
        medians[mode] = {
            k: float(np.median([s["timing"][k] for s in cost_rows if s["mode"] == mode]))
            for k in [
                "integration_quality_wall_s",
                "resident_body_s",
                "observer_output_s",
                "sync_readback_wall_s",
                "hdf5_output_s",
                "graph_capture_s",
                "cuda_event_trajectory_span_s",
            ]
        }
    processes = [dict(label=p.name, **js(p)) for p in sorted(BASE.glob("*_PROCESS.json"))]
    write(
        out / "COST_AND_YIELD.json",
        dict(
            matched_diagnostic_medians=medians,
            matched_rows=cost_rows,
            actual_science=trials + b2_trials,
            processes=processes,
            total_process_wall_s=sum(s["full_process_wall_s"] for s in processes),
            complete_10ms_samples=0,
            new_independent_initial_states=0,
            original_initial_states_reused=2,
            full_OFF_S_O_triplets=0,
            production_samples_per_hour=None,
            scope="Diagnostic costs; CUDA activity span includes host gaps. Small-window completions are not scientific production yield.",
        ),
    )
    gate = dict(
        status="ENGINEERING_COMPLETE_EXPLORATORY",
        LUB4_TWO_SHORT_WINDOWS_PASS=all(
            r["status"] == "COMPLETE" for r in trials if r["mode"].endswith("DT10") and r["target_s"] < 0.01
        ),
        LUB4_FULL_10MS_PASS=any(r["status"] == "COMPLETE" for r in trials if r["target_s"] == 0.01),
        LUBRICATION_NOT_FULLY_CONVERGED=True,
        MORPHOLOGY_NOT_CONVERGED=True,
        LUB4_TIMESTEP_CONVERGENCE="FAIL"
        if comps and any(c["status"] != "COMMON_PREFIX_PASS" for c in comps)
        else "NOT_EVALUATED",
        B2_GEOMETRIC_UNILATERALITY_PASS=bool(b2_rows)
        and all(not r["false_reactions"] and r["max_spring_law_force_error_N"] < 1e-15 for r in b2_rows),
        B2_TIMESTEP_CONVERGENCE="FAIL"
        if b2_comps and any(c["status"] != "COMMON_PREFIX_PASS" for c in b2_comps)
        else "B2_TIMESTEP_CONVERGENCE_NOT_EVALUATED",
        PRODUCTION_NOT_ADMITTED=True,
        PHYSICS_UNCALIBRATED=True,
        training_ready=False,
        default_unchanged=True,
        source=identity(),
    )
    if b2_comps and all(c["status"] == "COMMON_PREFIX_PASS" for c in b2_comps):
        gate["B2_TIMESTEP_CONVERGENCE"] = "COMMON_PREFIX_PASS_ONLY"
    write(out / "FINAL_GATE.json", gate)
    print("B2", [(c["maxima"], c["gates"]) for c in b2_comps])
    print([(r["case"], r["mode"], r["status"], r["accepted_steps"]) for r in trials])
    print([(c["case"], c["maxima"], c["gates"]) for c in comps])


if __name__ == "__main__":
    main()

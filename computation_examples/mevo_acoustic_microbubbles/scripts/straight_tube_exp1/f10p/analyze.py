"""Independent BASIC, COM, shape, topology and fixed-ID numerical admissions."""

import numpy as np
import h5py
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scripts.straight_tube_exp1.f10o.analyze import load, metrics, topology_record, compare
from .common import BASE, OLD, ROOT, CONFIG, read, write


def comparison(a, b, kind):
    summary, arrays = compare(a, b, kind)
    ta, xa = load(a)
    tb, xb = load(b)
    pairs = [(i, int(np.argmin(abs(tb - t)))) for i, t in enumerate(ta) if np.min(abs(tb - t)) < 1e-13]
    ma = metrics(ta[[i for i, j in pairs]], xa[[i for i, j in pairs]])
    mb = metrics(tb[[j for i, j in pairs]], xb[[j for i, j in pairs]])
    Qscale = np.linalg.norm(ma["tensor_m2"], axis=(-2, -1))
    arrays["Q_relative"] = arrays["tensor_m2"] / Qscale
    arrays["aspect_relative"] = arrays["aspect"] / ma["aspect"]
    distribution = []
    for x, y in zip(ma["centered_positions_m"], mb["centered_positions_m"]):
        distance = cdist(x, y, metric="sqeuclidean")
        i, j = linear_sum_assignment(distance)
        distribution.append(np.sqrt(distance[i, j].mean()) / read(CONFIG)["radius_m"])
    arrays["matched_RMS_over_R"] = np.array(distribution)
    identifiable = ma["axis_identifiable"] & mb["axis_identifiable"]
    angles = np.degrees(
        np.arccos(np.clip(abs(np.einsum("ti,ti->t", ma["principal_axis"], mb["principal_axis"])), 0, 1))
    )
    arrays["axis_angle_deg"] = np.where(identifiable, angles, np.nan)
    arrays["axis_identifiable_both"] = identifiable
    summary.update(
        {k: float(arrays[k].max(initial=0)) for k in ["Q_relative", "aspect_relative", "matched_RMS_over_R"]}
    )
    summary["axis_angle_deg"] = float(angles[identifiable].max(initial=0)) if identifiable.any() else None
    summary["axis_direction_scope"] = (
        "Only frames with both relative top eigenvalue gaps >0.001; degenerate directions NOT_IDENTIFIABLE"
    )
    return summary, arrays


def layered(record, dt, cache, topologies, policy):
    complete = record["status"] == "PASS" and abs(record["completed_s"] - 0.002) < 1e-13
    available = dt is not None and cache is not None
    comparisons = [dt, cache] if available else []
    full = available and all(c["complete_2ms"] for c in comparisons)
    position = complete and full and all(c["COM_m"] <= policy["group_COM_tolerance_m"] for c in comparisons)
    signal = float(np.linalg.norm(record["displacement_m"]))
    noise = max([c["COM_m"] for c in comparisons], default=float("inf"))
    signal_verified = signal > noise
    # OFF is a numerically admitted null reference, not demonstrated actuation.
    position = position and (signal_verified or record["source"] == "OFF")
    shape = (
        complete
        and full
        and all(
            c["Rg_m"] <= policy["group_Rg_tolerance_m"]
            and c["directional_xy"] <= policy["group_xy_tolerance"]
            and c["Q_relative"] <= policy["group_Q_relative"]
            and c["aspect_relative"] <= policy["group_aspect_relative"]
            and c["matched_RMS_over_R"] <= policy["distribution_RMS_over_R"]
            and (c["axis_angle_deg"] is None or c["axis_angle_deg"] <= policy["principal_axis_angle_deg"])
            for c in comparisons
        )
    )
    particle = (
        complete and full and all(c["fixed_id_m"] <= policy["fixed_id_tolerance_m"] for c in comparisons)
    )
    # Conservative gate: this initial implementation admits survival of the initial
    # single network only. Split/merge cases retain ledgers but require explicit review.
    topo = (
        complete
        and full
        and len(topologies) == 3
        and all(t["initial_network_persists_entire_2ms"] for t in topologies)
    )
    return dict(
        BASIC_NUMERICS_PASS="PASS" if complete else "FAIL",
        GROUP_POSITION_NUMERICALLY_ADMITTED="PASS" if position else ("FAIL" if available else "NOT_RUN"),
        GROUP_SHAPE_NUMERICALLY_ADMITTED="PASS" if shape else ("FAIL" if available else "NOT_RUN"),
        GROUP_TOPOLOGY_NUMERICALLY_ADMITTED="PASS" if topo else ("FAIL" if available else "NOT_RUN"),
        PARTICLE_ID_NUMERICALLY_ADMITTED="PASS" if particle else ("FAIL" if available else "NOT_RUN"),
        particle_status="PARTICLE_ID_NUMERICALLY_ADMITTED" if particle else "PARTICLE_ID_NOT_ADMITTED",
        simulation_group_position_eligible=bool(position),
        simulation_group_shape_eligible=bool(shape),
        simulation_topology_eligible=bool(topo),
        simulation_particle_id_eligible=bool(particle),
        control_signal_m=signal,
        COM_noise_m=noise if np.isfinite(noise) else None,
        control_signal_exceeds_numerical_difference=bool(signal_verified),
        topology_scope="Initial 32-member network survival for complete 2ms; later events may be right-censored. No physical fusion.",
        EXPERIMENTALLY_CALIBRATED="NOT_VERIFIED",
        TRAINING_READY=False,
        training_ready=False,
    )


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    p = read(CONFIG)
    runs = {}
    paths = {}
    comps = {}
    costs = []
    for campfile in sorted(BASE.glob("*/CAMPAIGN.json")):
        c = read(campfile)
        costs.append(c)
        for f in sorted((campfile.parent / "runs").glob("*/RESULT.json")):
            name = (
                campfile.parent.name + "__" if campfile.parent.name.startswith("pilot") else ""
            ) + f.parent.name
            r = read(f)
            paths[name] = f.parent
            t, x = load(f.parent)
            m = metrics(t, x)
            topo, labels, ets, eids = topology_record(t, x)
            np.savez_compressed(out / (name + ".npz"), **m, labels=labels, edge_times_s=ets, edge_ids=eids)
            write(out / (name + "_TOPOLOGY.json"), topo)
            with h5py.File(f.parent / "trajectory.h5") as h:
                rows = h["step_metrics"][:]
                pressure = abs(h["force_samples/pressure_Pa"][:])
            r.update(
                rg_initial_final_m=m["rg_m"][[0, -1]],
                xy_initial_final=m["directional_xy"][[0, -1]],
                aspect_initial_final=m["aspect"][[0, -1]],
                tensor_initial_final_m2=m["tensor_m2"][[0, -1]],
                angle_to_x_initial_final_deg=m["angle_to_x_deg"][[0, -1]],
                topology=topo,
                CPU_solver_s=float(rows[:, 7].sum()),
                CPU_resistance_s=float(rows[:, 8].sum()),
                measured_step_s=float(rows[:, 6].sum()),
                pressure_min_max_Pa=[float(pressure.min()), float(pressure.max())],
            )
            runs[name] = r
    groups = ["STRONG", "OFF", "WEAK"] + [
        name[:-4] for name in runs if name.startswith("NAV_") and name.endswith("_DT1")
    ]
    admissions = {}
    for tag in groups:
        name = tag + "_DT1"
        if name not in paths:
            continue
        topo = [runs[name]["topology"]]
        for suffix in ["DT05", "FINE"]:
            key = tag + "_" + suffix
            if key in paths:
                summary, arrays = comparison(paths[name], paths[key], "Same hard-sphere model " + suffix)
                comps[key] = summary
                np.savez_compressed(out / ("compare_" + key + ".npz"), **arrays)
                topo.append(runs[key]["topology"])
        admissions[tag] = layered(runs[name], comps.get(tag + "_DT05"), comps.get(tag + "_FINE"), topo, p)
    for tag, old in [("STRONG", "OPTIMIZED_REPLAY"), ("WEAK", "SHAPE_REPLAY")]:
        if tag + "_DT1" not in paths:
            continue
        for model, folder in [
            ("B2", ROOT / "results/straight_tube_exp1/f10o/validation_v1/runs" / (tag + "_DT1")),
            ("M0", OLD / "campaign_v1/runs" / old),
        ]:
            key = tag + "_" + model + "_VS_HARD"
            s, a = comparison(folder, paths[tag + "_DT1"], "Different contact models; not integration error")
            comps[key] = s
            np.savez_compressed(out / ("compare_" + key + ".npz"), **a)
    for suffix in ["DT05", "FINE"]:
        name = "STRONG_PREFIX_" + suffix
        if name in paths and "STRONG_DT1" in paths:
            s, a = comparison(
                paths["STRONG_DT1"], paths[name], "Same model; common legal prefix only, NOT 2ms admission"
            )
            comps[name] = s
            np.savez_compressed(out / ("compare_" + name + ".npz"), **a)
    # Instantaneous edge timing remains separate from initial network survival.
    edge_checks = {}
    for tag in ["OFF", "WEAK"]:
        for suffix in ["DT05", "FINE"]:
            if tag + "_DT1" not in paths or tag + "_" + suffix not in paths:
                continue
            ledgers = []
            for name in [tag + "_DT1", tag + "_" + suffix]:
                d = np.load(out / (name + ".npz"))
                ledger = {}
                for t, e in zip(d["edge_times_s"], d["edge_ids"]):
                    ledger.setdefault(tuple(e.tolist()), []).append(float(t))
                ledgers.append(ledger)
            a, b = ledgers
            matched = a.keys() == b.keys() and all(len(a[k]) == len(b[k]) for k in a)
            error = (
                max([max(abs(np.array(a[k]) - b[k]), default=0) for k in a], default=0) if matched else None
            )
            edge_checks[tag + "_" + suffix] = dict(
                matched_event_sequences=matched,
                max_time_error_s=error,
                status="PASS" if matched and error <= p["instant_edge_time_tolerance_s"] else "FAIL",
                event_counts=[sum(map(len, a.values())), sum(map(len, b.values()))],
                scope="Instantaneous 2.1R edges, not persistence/fusion",
            )
    gate = dict(
        labels=p["labels"],
        ENGINEERING_PASS="PENDING_FINAL_CHECKS",
        groups=admissions,
        navigation_admitted=bool(
            admissions.get("STRONG", {}).get("simulation_group_position_eligible", False)
        ),
        targets_m=p["XY_targets_m"],
        target_tolerance_m=p["XY_tolerance_m"],
        TARGET_REACHED={},
        EXPERIMENTALLY_CALIBRATED="NOT_VERIFIED",
        TRAINING_READY=False,
        training_ready=False,
        historical_B2_strong="Rejected at 0.175ms; never extrapolated",
        windows_5_10_30ms="NOT_RUN",
    )
    for tag, g in admissions.items():
        r = runs[tag + "_DT1"]
        end = np.asarray(r["displacement_m"])
        g["target_error_m"] = {
            axis: float(np.linalg.norm(end - np.asarray(target)))
            for axis, target in p["XY_targets_m"].items()
        }
        g["TARGET_REACHED"] = {
            axis: bool(g["simulation_group_position_eligible"] and e <= p["XY_tolerance_m"])
            for axis, e in g["target_error_m"].items()
        }
        g["SHAPE_TARGET_REACHED"] = bool(
            g["simulation_group_shape_eligible"] and r["xy_initial_final"][-1] >= p["shape_target_xy"]
        )
        gate["TARGET_REACHED"][tag] = g["TARGET_REACHED"]
    gate["group_pass_particle_fail_actual_cases"] = [
        tag
        for tag, g in admissions.items()
        if g["simulation_group_position_eligible"] and not g["simulation_particle_id_eligible"]
    ]
    local = (
        read(BASE / "local_foundations_v1.json")
        if (BASE / "local_foundations_v1.json").exists()
        else {"wall_s": 0}
    )
    local["wall_s"] += read(BASE / "frozen_failure_recheck_v2.json")["wall_s"]
    gate["instantaneous_edge_checks"] = edge_checks
    old_dt = read(ROOT / "results/straight_tube_exp1/f10o/analysis/COMPARISONS.json")["WEAK_DT05"]
    gap = comps["WEAK_B2_VS_HARD"]["Rg_m"]
    noise = max(comps["WEAK_DT05"]["Rg_m"], comps["WEAK_FINE"]["Rg_m"], old_dt["Rg_m"])
    gate["CONTACT_MODEL_DEPENDENT"] = bool(gap > 10 * noise)
    gate["model_Rg_difference_vs_max_numerical_difference"] = dict(model_m=gap, numerical_m=noise)
    gate["contact_model_scope"] = (
        "M0, B2 and hard-sphere weak-control shape changes differ; see paired actual metrics, not pure integration error."
    )
    gate["stop_reason"] = (
        "Strong local QP passes but actual nonlinear-gap complementarity fails at 0.464ms; no full 2ms strong admission or XY search."
    )
    gate["navigation_status"] = "NOT_RUN_STRONG_BASIC_QUALITY_FAILED"
    write(out / "RUN_ANALYSIS.json", runs)
    write(out / "COMPARISONS.json", comps)
    write(out / "FINAL_GATE.json", gate)
    write(
        out / "COSTS.json",
        dict(
            remote_process_s=sum(c["wall_s"] for c in costs),
            local_foundations_s=local["wall_s"],
            total_numerical_s=sum(c["wall_s"] for c in costs) + local["wall_s"],
            integration_s=sum(r["integration_s"] for r in runs.values()),
            CPU_solver_s=sum(r["CPU_solver_s"] for r in runs.values()),
            CPU_resistance_s=sum(r["CPU_resistance_s"] for r in runs.values()),
            GPU_kernel_only_s=None,
            GPU_note="GPU field/force/walls used; pure kernel time not separately profiled",
            peak_RSS_MiB=max([c.get("peak_rss_MiB", 0) for c in costs], default=0),
        ),
    )
    print(
        {
            k: {s: v for s, v in g.items() if s.startswith(("GROUP", "PARTICLE", "BASIC"))}
            for k, g in admissions.items()
        }
    )


if __name__ == "__main__":
    main()

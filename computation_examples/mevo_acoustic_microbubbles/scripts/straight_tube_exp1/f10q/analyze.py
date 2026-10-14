"""Reuse F10P layered admission; accepted prefixes never extended past rejection."""

import numpy as np
import h5py
from scripts.straight_tube_exp1.f10o.analyze import load, metrics, topology_record
from scripts.straight_tube_exp1.f10p.analyze import comparison, layered
from .common import BASE, ROOT, CONFIG, read, write, identity


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    p = read(CONFIG)
    runs = {}
    paths = {}
    costs = []
    comps = {}
    admissions = {}
    for cp in sorted(BASE.glob("*/CAMPAIGN.json")):
        camp = read(cp)
        costs.append(camp)
        for rp in sorted((cp.parent / "runs").glob("*/RESULT.json")):
            tag = rp.parent.name
            r = read(rp)
            paths[tag] = rp.parent
            t, x = load(rp.parent)
            m = metrics(t, x)
            top, labs, ets, eids = topology_record(t, x)
            np.savez_compressed(out / (tag + ".npz"), **m, labels=labs, edge_times_s=ets, edge_ids=eids)
            with h5py.File(rp.parent / "trajectory.h5") as h:
                rows = h["step_metrics"][:]
                extra = h["extra_metrics"][:]
                mask = rows[:, 0].astype(bool)
                accepted = rows[mask]
                # All state nodes must also be legal in the true central cavity branch.
                radial = np.linalg.norm(x[:, :, 1:], axis=2)
                assert (
                    radial.max() + p["radius_m"] < 0.0008 and abs(x[:, :, 0]).max() + p["radius_m"] < 0.01625
                )
            r.update(
                topology=top,
                rg_initial_final_m=m["rg_m"][[0, -1]],
                xy_initial_final=m["directional_xy"][[0, -1]],
                aspect_initial_final=m["aspect"][[0, -1]],
                Q_initial_final_m2=m["tensor_m2"][[0, -1]],
                direction_initial_final_deg=m["angle_to_x_deg"][[0, -1]],
                CPU_solver_s=float(rows[:, 7].sum()),
                CPU_resistance_s=float(rows[:, 8].sum()),
                pressure_min_max_Pa=[float(extra[:, 3].min()), float(extra[:, 4].max())],
                force_max_N=float(extra[:, 5].max()),
                accepted_max_overlap_m=float(accepted[:, 2].max(initial=0)),
                accepted_max_true_comp_m=float(accepted[:, 17].max(initial=0)),
                accepted_max_residual_m=float(accepted[:, 1].max(initial=0)),
                accepted_max_force_residual_N=float(accepted[:, 14].max(initial=0)),
                accepted_max_local_residual_m=float(extra[mask, 1].max(initial=0)),
                QP_count_distribution={
                    str(int(k)): int((rows[:, 9] == k).sum()) for k in np.unique(rows[:, 9])
                },
            )
            runs[tag] = r
    for tag in ["OFF", "STRONG", "NAV_X", "NAV_Y"]:
        key = tag + "_DT1"
        if key not in paths:
            continue
        tops = [runs[key]["topology"]]
        for suffix in ["DT05", "FINE"]:
            other = tag + "_" + suffix
            if other in paths:
                s, a = comparison(paths[key], paths[other], "Same F10Q SQP: " + suffix)
                comps[other] = s
                np.savez_compressed(out / ("compare_" + other + ".npz"), **a)
                tops.append(runs[other]["topology"])
        gate = layered(runs[key], comps.get(tag + "_DT05"), comps.get(tag + "_FINE"), tops, p)
        gate["CONTACT_MODEL_DEPENDENT"] = (
            "NOT_VERIFIED_FOR_THIS_NEW_CONTROL"
            if tag.startswith("NAV")
            else "KNOWN_FROM_HISTORICAL_CONTRASTS"
        )
        gate["target_error_m"] = {
            a: float(np.linalg.norm(np.array(runs[key]["displacement_m"]) - target))
            for a, target in p["XY_targets_m"].items()
        }
        gate["TARGET_REACHED"] = {
            a: bool(gate["simulation_group_position_eligible"] and error <= p["XY_tolerance_m"])
            for a, error in gate["target_error_m"].items()
        }
        gate["relative_Rg_change"] = float(
            runs[key]["rg_initial_final_m"][-1] / runs[key]["rg_initial_final_m"][0] - 1
        )
        gate["navigation_shape_cost_acceptable"] = bool(
            abs(gate["relative_Rg_change"]) <= 0.05 and runs[key]["topology"]["max_components"] == 1
        )
        admissions[tag] = gate
    if "STRONG_DT1" in paths:
        for name, folder in [
            ("OLD_HARD", ROOT / "results/straight_tube_exp1/f10p/validation_v1/runs/STRONG_DT1"),
            ("B2", ROOT / "results/straight_tube_exp1/f10o/validation_v1/runs/STRONG_DT1"),
        ]:
            s, a = comparison(
                folder, paths["STRONG_DT1"], "Different contact discretizations on exact common legal times"
            )
            comps[name + "_VS_SQP"] = s
            np.savez_compressed(out / ("compare_" + name + "_VS_SQP.npz"), **a)
    # Full scientific budgets include local finite numerical work, not just remote integration.
    localA = sum(
        read(BASE / rel)["wall_s"]
        for rel in [
            "analysis/TRUE_GAP_FAILURE_DIAGNOSIS.json",
            "local/FOUNDATIONS_v1.json",
            "local/FOUNDATIONS_v2.json",
        ]
    )
    localB = sum(
        read(BASE / rel)["wall_s"] for rel in ["analysis/THREE_MODEL_COMPARISON.json", "CANDIDATES.json"]
    )
    route = {a: sum(c["wall_s"] for c in costs if c["route"] == a) for a in ["A", "B", "reserve"]}
    route["A"] += localA
    route["B"] += localB
    totals = dict(
        routes_wall_s=route,
        total_numerical_s=sum(route.values()),
        remote_process_s=sum(c["wall_s"] for c in costs),
        local_A_s=localA,
        local_B_s=localB,
        CPU_solver_s=sum(r["CPU_solver_s"] for r in runs.values()),
        CPU_resistance_s=sum(r["CPU_resistance_s"] for r in runs.values()),
        integration_s=sum(r["integration_s"] for r in runs.values()),
        peak_RSS_MiB=max(c.get("peak_rss_MiB", 0) for c in costs),
        GPU_kernel_only_s=None,
        GPU_scope="Existing Warp forces/field/walls used; CPU solver plus transfers and Python overhead. No pure GPU kernel timing claim.",
    )
    gate = dict(
        identity=identity(),
        labels=p["labels"],
        ENGINEERING_PASS="PENDING_FINAL_CHECKS",
        routes_independent=True,
        groups=admissions,
        new_candidates_evaluated=[k for k in runs if "_P" in k],
        short_candidate_scope="PASS is own 0.5ms basic quality only; not full 2ms layered admission",
        CONTACT_MODEL_DEPENDENT=True,
        contact_scope="Historical common-time shape contrasts establish model dependence; new navigation has only same-model numerical checks.",
        EXPERIMENTALLY_CALIBRATED="NOT_VERIFIED",
        TRAINING_READY=False,
        training_ready=False,
        windows_5_10_30ms="NOT_RUN",
        group_pass_particle_fail_actual_cases=[
            k
            for k, v in admissions.items()
            if v["simulation_group_position_eligible"] and not v["simulation_particle_id_eligible"]
        ],
        position_pass_shape_fail_actual_cases=[
            k
            for k, v in admissions.items()
            if v["simulation_group_position_eligible"] and not v["simulation_group_shape_eligible"]
        ],
    )
    gate["candidate_admissions"] = {}
    for name, record in runs.items():
        if "_P" not in name:
            continue
        gate["candidate_admissions"][name] = dict(
            BASIC_NUMERICS_PASS="PASS" if record["status"] == "PASS" else "FAIL",
            window_s=record["requested_s"],
            GROUP_POSITION_NUMERICALLY_ADMITTED="NOT_RUN",
            GROUP_SHAPE_NUMERICALLY_ADMITTED="NOT_RUN",
            GROUP_TOPOLOGY_NUMERICALLY_ADMITTED="NOT_RUN",
            PARTICLE_ID_NUMERICALLY_ADMITTED="NOT_RUN",
            TARGET_REACHED=False,
            CONTACT_MODEL_DEPENDENT="NOT_VERIFIED_FOR_THIS_CONTROL",
            EXPERIMENTALLY_CALIBRATED="NOT_VERIFIED",
            TRAINING_READY=False,
            scope="Short-window own quality only; no dt/cache repeats of this short record",
        )
    write(out / "RUN_ANALYSIS.json", runs)
    write(out / "COMPARISONS.json", comps)
    write(out / "FINAL_GATE.json", gate)
    write(out / "COSTS.json", totals)
    print("costs", totals)
    print("admissions", admissions)


if __name__ == "__main__":
    main()

"""Postprocess completed real runs, preserving all incomplete and failed outcomes."""

from collections import Counter
from pathlib import Path
import shutil
import h5py
import numpy as np
from scipy.spatial.distance import pdist
from scripts.stage5a_r3.common import BASE, R2, read, write, sha
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r3.lineage import compare_all
from scripts.stage5a_r3.current import campaign as current_campaign


def cluster_report(comparison, original_group):
    arrays = np.load(comparison["samples"])
    result = {}
    for label, mask in [
        ("all_active", np.ones(len(arrays["particle_id"]), bool)),
        ("original_group", np.isin(arrays["particle_id"], original_group)),
        ("complement", ~np.isin(arrays["particle_id"], original_group)),
    ]:
        a, b = [arrays[k][:, mask] for k in ["left_positions_m", "right_positions_m"]]
        ac, bc = a.mean(1), b.mean(1)
        ar = np.sqrt(np.mean(np.sum((a - ac[:, None]) ** 2, axis=2), axis=1))
        br = np.sqrt(np.mean(np.sum((b - bc[:, None]) ** 2, axis=2), axis=1))
        aq, bq = [np.quantile(pdist(x[-1]), [0.25, 0.5, 0.75]) for x in [a, b]]
        travel = float(np.linalg.norm(np.diff(ac, axis=0), axis=1).sum())
        budget = max(0.1 * float(arrays["radii_m"][0]), 0.02 * travel)
        error = float(np.linalg.norm(bc[-1] - ac[-1]))
        rg_end = float(abs(br[-1] / ar[-1] - 1))
        rg_series = float(np.linalg.norm(br - ar) / np.linalg.norm(ar))
        quart = np.abs(bq / aq - 1)
        bound = arrays["left_sampling_bound_m"] + arrays["right_sampling_bound_m"]
        uncertainties = np.r_[
            bound[-1] / ar[-1], np.linalg.norm(bound) / np.linalg.norm(ar), 2 * bound[-1] / aq, bound[-1]
        ]
        observed = np.r_[rg_end, rg_series, quart, error]
        limits = np.r_[np.full(5, 0.02), budget]
        # A conservative unknown sampling interval never silently passes a boundary.
        status = (
            "NOT_EVALUABLE"
            if not comparison.get("sampling_bound_certified", False)
            else "PASS"
            if np.all(observed + uncertainties <= limits)
            else "FAIL"
            if np.any(observed - uncertainties > limits)
            else "NOT_EVALUABLE"
        )
        result[label] = dict(
            status=status,
            Rg_endpoint_relative=rg_end,
            Rg_series_relative_L2=rg_series,
            pair_quartiles_relative=quart.tolist(),
            centroid_difference_m=error,
            centroid_budget_m=budget,
            reference_sampled_centroid_travel_m=travel,
            endpoint_sampling_bound_m=float(bound[-1]),
            scope="Descriptive consistency against the named numerical reference, not independent physical accuracy",
        )
    return result


def add_comparison(case, name, left, right, group):
    left, right = local(left), local(right)
    if not left.exists() or not right.exists():
        return dict(
            case=case,
            name=name,
            status="NOT_RUN",
            reason="LOCAL_SOURCE_NOT_SYNCED",
            left=str(left),
            right=str(right),
        )
    cache = BASE / "comparisons" / f"{case}_{name}_REPORT.json"
    hashes = [sha(left), sha(right)]
    if cache.exists():
        previous = read(cache)
        if (
            [previous["left_sha256"], previous["right_sha256"]] == hashes
            and previous.get("sampling_audit_version") == 2
            and previous.get("neighbor_trace_version") == 2
        ):
            return previous
    result = compare_all(case, name, left, right, group)
    result["reason"] = (
        "No independent all-active global error bound established; descriptive fixed-ID comparison only"
    )
    result["cluster_consistency"] = cluster_report(result, group)
    arrays = np.load(result["samples"])
    errors = arrays["all_active_raw_errors_m"]
    j, i = np.unravel_index(errors.argmax(), errors.shape)
    ids = arrays["particle_id"]
    # Save the same identified bubble and true nearest neighbors in each method.
    a, b = arrays["left_positions_m"], arrays["right_positions_m"]
    near = np.unique(
        np.r_[
            np.argsort(np.linalg.norm(a[j] - a[j, i], axis=1))[:8],
            np.argsort(np.linalg.norm(b[j] - b[j, i], axis=1))[:8],
        ]
    )
    output = BASE / "comparisons" / f"{case}_{name}_WORST_NEIGHBORS.npz"
    full_traces = {}
    for label, path in [("left", left), ("right", right)]:
        with h5py.File(path) as h:
            index = np.array([int(np.flatnonzero(h["particle_id"][:] == pid)[0]) for pid in ids[near]])
            full_traces[label + "_all_accepted_time_s"] = h["time_s"][:]
            full_traces[label + "_all_accepted_positions_m"] = h["positions_m"][:, index]
    np.savez_compressed(
        output,
        particle_id=ids[near],
        worst_particle_id=ids[i],
        clock_s=arrays["clock_s"],
        left_actual_time_s=arrays["left_actual_time_s"],
        right_actual_time_s=arrays["right_actual_time_s"],
        left_positions_m=a[:, near],
        right_positions_m=b[:, near],
        radii_m=arrays["radii_m"][near],
        **full_traces,
    )
    result["worst_neighborhood"] = str(output)
    result["neighbor_trace_version"] = 2
    result["worst_neighborhood_sha256"] = sha(output)
    radius = float(arrays["radii_m"][0])
    targets = dict(RMS_m=0.01 * radius, max_m=0.05 * radius)
    series_rms = np.sqrt(np.mean(errors**2, axis=1))
    result["self_resolution_targets"] = targets
    result["resolution_comparison"] = dict(
        status="PASS"
        if series_rms.max() <= targets["RMS_m"] and errors.max() <= targets["max_m"]
        else "FAIL",
        max_series_RMS_m=float(series_rms.max()),
        max_error_m=float(errors.max()),
        scope="Predeclared numerical resolution targets only; accuracy status remains NOT_EVALUABLE",
    )
    write(cache, result)
    return result


def method_usage(folder):
    accepted = read(folder / "ACCEPTED.jsonl")
    rejected = read(folder / "REJECTED.jsonl")
    summary = read(folder / "SUMMARY.json")
    with h5py.File(folder / "trajectories.h5") as h:
        times = h["time_s"][:]
        masks = h["active"][:]
        ids, radii = h["particle_id"][:], h["radii_m"][:]
        finite_active = all(np.isfinite(h["positions_m"][i][mask]).all() for i, mask in enumerate(masks))
        actual_clock = bool(
            len(times) == len(accepted) + 1
            and np.all(np.diff(times) > 0)
            and np.array_equal(times[1:], [r["time_s"] for r in accepted])
            and abs(sum(r["dt_s"] for r in accepted) - (times[-1] - times[0])) < 2e-12
        )
        checks = dict(
            finite_active_positions=bool(finite_active),
            strictly_increasing_actual_time=bool(np.all(np.diff(times) > 0)),
            accepted_time_ledger_matches=actual_clock,
            fixed_active_mask=bool(np.all(masks == masks[0])),
            unique_particle_ids=bool(len(np.unique(ids)) == len(ids)),
            positive_finite_radii=bool(np.isfinite(radii).all() and np.all(radii > 0)),
            full_expected_active_count=int(masks[0].sum()) == {"B2": 186, "B3": 189}[summary["case"]],
            nonempty_accepted_history=bool(accepted),
        )
        integrity = dict(
            status="PASS" if all(checks.values()) else "FAIL",
            checks=checks,
            active_count=int(masks[0].sum()),
            actual_start_s=float(times[0]),
            actual_end_s=float(times[-1]),
            trajectory_sha256=sha(folder / "trajectories.h5"),
            inactive_sentinels="Inactive NaN positions are not active-state failures",
        )
    counts = []
    for method in ["S0", "S1"]:
        rows = [r for r in accepted if r["substep_method"] == method]
        bursts = []
        current = []
        for r in accepted:
            if r["substep_method"] == method:
                current.append(r)
            elif current:
                bursts.append(current)
                current = []
        if current:
            bursts.append(current)
        counts.append(
            dict(
                method=method,
                accepted_substeps=len(rows),
                physical_duration_s=sum(r["dt_s"] for r in rows),
                dt_quantiles_s=np.quantile([r["dt_s"] for r in rows], [0, 0.5, 1]).tolist() if rows else None,
                bursts=[
                    dict(
                        start_s=x[0]["time_s"] - x[0]["dt_s"],
                        end_s=x[-1]["time_s"],
                        accepted_substeps=len(x),
                        duration_s=sum(r["dt_s"] for r in x),
                    )
                    for x in bursts
                ],
                nonlinear_kkt_max=max(
                    (r["nonlinear_kkt"] for r in rows if r["nonlinear_kkt"] is not None), default=None
                ),
                nonlinear_kkt_status="NOT_APPLICABLE" if method == "S0" else "EVALUATED",
                accepted_nonlinear_iterations=sum(r["nonlinear_iterations"] for r in rows),
                accepted_qp_iterations=sum(r["qp_iterations"] for r in rows),
                accepted_native_fallbacks=sum(r["fallbacks"] for r in rows),
            )
        )
    return dict(
        run=str(folder),
        case=summary["case"],
        target_s=summary["target_s"],
        status=summary["status"],
        methods=counts,
        state_integrity=integrity,
        continuous_certificate_all_pass=all(r["safety"]["pass_gate"] for r in accepted),
        accepted_count=len(accepted),
        rejected_count=len(rejected),
        rejection_reasons=dict(Counter(r["reason"] for r in rejected)),
        rejected_trial_wall_s=sum(r["trial_wall_s"] for r in rejected),
        accepted_safety=dict(
            kkt_max=summary["accepted_kkt_max"],
            nonlinear_kkt_max=summary["accepted_nonlinear_kkt_max"],
            linear_residual_max=summary["accepted_linear_residual_max"],
            pair_swept_gap_min_m=summary["accepted_pair_swept_gap_min_m"],
            wall_swept_bound_min_m=summary["accepted_wall_swept_gap_bound_min_m"],
            projection_max_m=summary["projection_max_m"],
        ),
        rejected_scope="Rejected candidates excluded from all accepted-state residual and safety maxima",
        failed_trial_candidates=[r for r in rejected if r["reason"] != "LOCAL_ERROR"],
        cost=dict(
            full_run_wall_s=summary["wall_s"],
            exclusive_s=summary["timings_exclusive_s"],
            qp_calls=summary["qp_calls"],
            IO_inclusive_s=summary["io_seconds"],
            monitor_inclusive_s=summary["monitor_inclusive_seconds"],
        ),
    )


def main():
    campaign = current_campaign()
    old = read(R2 / "final/CAMPAIGN_V4A/RESULT.json")["ledger"]
    lineage = read(BASE / "final/REFERENCE_LINEAGE_AUDIT.json")
    ledger = campaign["ledger"]
    comparisons = []
    usages = []
    for case in ["B2", "B3"]:
        group = np.load(R2 / "comparisons" / f"{case}_0.001_v4a_legacy_reference.npz")["particle_ids"]
        archived = next(r["paths"]["archived"] for r in lineage["rows"] if r["case"] == case)
        for duration in [0.001, 0.005]:
            entries = {
                r["kind"]: r
                for r in ledger
                if r["case"] == case and r["duration_s"] == duration and r["repeat"] == 0 and r.get("run")
            }
            if "nominal" not in entries:
                continue
            candidate = local(entries["nominal"]["run"]) / "trajectories.h5"
            prefix = f"{campaign.get('version', 'R3_V1')}_{duration:g}"
            for kind in ["refined", "reference"]:
                if kind in entries:
                    comparisons.append(
                        add_comparison(
                            case,
                            prefix + "_" + kind,
                            local(entries[kind]["run"]) / "trajectories.h5",
                            candidate,
                            group,
                        )
                    )
            if "reference" in entries:
                previous_fresh = next(
                    (
                        r
                        for r in old
                        if r["case"] == case
                        and r["duration"] == duration
                        and r["kind"] == "reference"
                        and r.get("repeat") == 0
                        and r.get("run")
                    ),
                    None,
                )
                if previous_fresh:
                    comparisons.append(
                        add_comparison(
                            case,
                            prefix + "_old_fresh_vs_new_fresh",
                            Path(previous_fresh["run"]) / "trajectories.h5",
                            local(entries["reference"]["run"]) / "trajectories.h5",
                            group,
                        )
                    )
            # Both legacy lineages retained. Archived 5 ms is compared only to 5 ms.
            if duration == 0.005:
                comparisons.append(add_comparison(case, prefix + "_archived", archived, candidate, group))
                if "reference" in entries:
                    comparisons.append(
                        add_comparison(
                            case,
                            prefix + "_archived_vs_fresh",
                            archived,
                            local(entries["reference"]["run"]) / "trajectories.h5",
                            group,
                        )
                    )
            if case == "B3" and campaign.get("version") == "R3_V2":
                previous = read(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")["ledger"]
                former = next(
                    r
                    for r in previous
                    if r["case"] == case and r["kind"] == "nominal" and r["duration_s"] == duration
                )
                comparisons.append(
                    add_comparison(
                        case,
                        prefix + "_V1_global_vs_V2_serial",
                        local(former["run"]) / "trajectories.h5",
                        candidate,
                        group,
                    )
                )
            oldrow = next(
                (
                    r
                    for r in old
                    if r["case"] == case
                    and r["duration"] == duration
                    and r["kind"] == "nominal"
                    and r.get("repeat") == 0
                    and r.get("run")
                ),
                None,
            )
            if oldrow:
                comparisons.append(
                    add_comparison(
                        case, prefix + "_R2", Path(oldrow["run"]) / "trajectories.h5", candidate, group
                    )
                )
            if duration == 0.005:
                for repeat in [1, 2]:
                    other = next(
                        (
                            r
                            for r in ledger
                            if r["case"] == case
                            and r["duration_s"] == duration
                            and r["repeat"] == repeat
                            and r["kind"] == "nominal"
                            and r.get("run")
                        ),
                        None,
                    )
                    if other:
                        comparisons.append(
                            add_comparison(
                                case,
                                prefix + f"_repeat{repeat}",
                                candidate,
                                local(other["run"]) / "trajectories.h5",
                                group,
                            )
                        )
    for entry in ledger:
        if entry.get("run") and entry["kind"] in ["nominal", "refined"]:
            folder = local(entry["run"])
            if (folder / "SUMMARY.json").exists():
                usages.append(method_usage(folder))
    prior = BASE / "final/ALL_ACTIVE_TRAJECTORY_AUDIT_PRE_R3.json"
    current = BASE / "final/ALL_ACTIVE_TRAJECTORY_AUDIT.json"
    if not prior.exists():
        shutil.copy2(current, prior)
    write(
        current,
        dict(
            status="NOT_EVALUABLE",
            reason="Independent all-active global bound not established; no automatic promotion",
            pre_R3=read(prior),
            formal=comparisons,
            early_B2_local_precision=read(BASE / "final/B2_EARLY_PRECISION_AUDIT.json")
            if (BASE / "final/B2_EARLY_PRECISION_AUDIT.json").exists()
            else None,
            campaign_complete=campaign["complete"],
        ),
    )
    write(BASE / "final/FORMAL_RUN_AUDIT.json", dict(campaign_complete=campaign["complete"], rows=usages))
    method_path = BASE / "final/S0_S1_USAGE_AND_REFINEMENT.json"
    methods = read(method_path)
    methods["formal_usage"] = usages
    write(method_path, methods)
    cost_path = BASE / "final/DYNAMIC_COST_PROFILE.json"
    cost = read(cost_path)
    cost["formal_runs"] = [
        dict(run=r["run"], cost=r["cost"], rejected_trial_wall_s=r["rejected_trial_wall_s"]) for r in usages
    ]
    instrumented = []
    for path in (BASE / "diagnosis").glob("DYNAMIC_INSTRUMENTED_*.json"):
        row = read(path)
        instrumented.append(
            {k: v for k, v in row.items() if k not in ["trials", "stages", "qp_calls", "native_candidates"]}
        )
        instrumented[-1].update(
            raw_path=str(path),
            raw_sha256=sha(path),
            qp_count=len(row["qp_calls"]),
            native_candidate_count=len(row["native_candidates"]),
            mean_QP_adapter_wall_s=float(np.mean([q["wall_s"] for q in row["qp_calls"]]))
            if row["qp_calls"]
            else None,
            QP_adapter_wall_sum_s=sum(q["wall_s"] for q in row["qp_calls"]),
            native_candidate_wall_sum_s=sum(q["wall_s"] for q in row["native_candidates"]),
        )
    cost["instrumented_windows"] = instrumented
    write(cost_path, cost)
    write(
        BASE / "final/FORMAL_COMPARISONS.json", dict(campaign_complete=campaign["complete"], rows=comparisons)
    )
    lineage["formal_comparisons"] = comparisons
    lineage["formal_campaign_complete"] = campaign["complete"]
    lineage["historical_sampling_bound_clarification"] = (
        "Pre-R3 reports retain their original values and identities. Their neighboring-frame displacement labeled sampling bound is not a proven bound for unmatched legacy times with projection; treat it as a displacement scale. V2 postprocessing records a certification flag and withholds consistency gates when that guarantee is absent. Exact common-frame comparisons are unaffected."
    )
    write(BASE / "final/REFERENCE_LINEAGE_AUDIT.json", lineage)
    print(
        dict(comparisons=len(comparisons), audited_runs=len(usages), complete=campaign["complete"]),
        flush=True,
    )


if __name__ == "__main__":
    main()

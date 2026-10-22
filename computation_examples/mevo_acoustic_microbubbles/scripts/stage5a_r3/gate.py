"""Separate motion, accuracy, reference and measured performance verdicts."""

import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, read, write, sha
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r3.campaign import numerical
from scripts.stage5a_r3.current import campaign as current_campaign, protocol_path


def same_initial(left, right):
    paths = [local(r["run"]) / "trajectories.h5" for r in [left, right]]
    if not all(p.exists() for p in paths):
        return None
    with h5py.File(paths[0]) as a, h5py.File(paths[1]) as b:
        return bool(
            all(
                np.array_equal(
                    a[k][0] if k in ["time_s", "positions_m", "active"] else a[k][:],
                    b[k][0] if k in ["time_s", "positions_m", "active"] else b[k][:],
                    equal_nan=True,
                )
                for k in ["time_s", "positions_m", "active", "particle_id", "radii_m"]
            )
        )


def completed(row, duration):
    s = row.get("summary") or {}
    return s.get("status") == "COMPLETE" and abs(s.get("completed_duration_s", -1) - duration) < 1e-13


def comparison_uses_run(comparison, entry):
    """Do not attach repeat-0 precision evidence to a different executable run."""
    if not entry.get("run"):
        return False
    folder = local(entry["run"])
    return any(local(comparison[key]).parent == folder for key in ["left", "right"] if comparison.get(key))


def ratio(reference, candidate, duration):
    if (
        not reference
        or not candidate
        or not completed(reference, duration)
        or not completed(candidate, duration)
    ):
        return dict(
            status="NOT_EVALUABLE",
            speedup=None,
            reason="Both same-duration complete tasks are required; no extrapolation",
        )
    same = same_initial(reference, candidate)
    if same is not True:
        return dict(
            status="NOT_EVALUABLE",
            speedup=None,
            reason="Initial-state identity not verified",
            same_initial=same,
        )
    old = reference.get("charged_wall_s", reference.get("wall_s_with_process_start"))
    new = candidate.get("charged_wall_s", candidate.get("wall_s_with_process_start"))
    if old is None or new is None:
        return dict(
            status="NOT_EVALUABLE", speedup=None, reason="Comparable process-inclusive wall times missing"
        )
    value = old / new
    return dict(
        status="PASS" if value >= 5 else "FAIL",
        non_regression_status="PASS" if value >= 1 else "FAIL",
        speedup=value,
        reference_wall_s=old,
        candidate_wall_s=new,
        reference_run=reference["run"],
        candidate_run=candidate["run"],
        target=5,
        same_initial=True,
        scope="Completion wall time including process start, initialization, solve, scientific output and monitoring",
    )


def main():
    campaign = current_campaign()
    version = campaign.get("version", "R3_V1")
    ledger = campaign["ledger"]
    comparisons = read(BASE / "final/FORMAL_COMPARISONS.json")["rows"]
    old = read(R2 / "final/CAMPAIGN_V4A/RESULT.json")["ledger"]
    tests = read(BASE / f"final/TEST_GATE_{version}.json")
    holdout = read(BASE / f"final/HOLDOUT_{version}.json")
    audits = {str(local(r["run"])): r for r in read(BASE / "final/FORMAL_RUN_AUDIT.json")["rows"]}
    rows = []
    for entry in ledger:
        case, duration, kind, repeat = [entry[k] for k in ["case", "duration_s", "kind", "repeat"]]
        summary = entry.get("summary")
        current = entry["status"]
        numerical_status = (
            "PASS"
            if kind != "reference" and numerical(summary)
            else "NOT_EVALUABLE"
            if kind == "reference"
            else "NOT_RUN"
            if current == "NOT_RUN"
            else "INCOMPLETE"
            if current in ["RUNNING", "INCOMPLETE"]
            else "FAIL"
        )
        state_audit = audits.get(str(local(entry["run"]))) if entry.get("run") else None
        if numerical_status == "PASS" and kind in ["nominal", "refined"]:
            if not state_audit or "state_integrity" not in state_audit:
                numerical_status = "NOT_EVALUABLE"
            elif (
                state_audit["state_integrity"]["status"] != "PASS"
                or not state_audit["continuous_certificate_all_pass"]
            ):
                numerical_status = "FAIL"
        comp = [
            r
            for r in comparisons
            if r["case"] == case
            and r["name"].startswith(f"{version}_{duration:g}_")
            and r.get("cluster_consistency")
            and comparison_uses_run(r, entry)
        ]
        ref = next(
            (
                r
                for r in ledger
                if r["case"] == case
                and r["duration_s"] == duration
                and r["kind"] == "reference"
                and r["repeat"] == repeat
            ),
            None,
        )
        protective = next(
            (
                r
                for r in ledger
                if r["case"] == case
                and r["duration_s"] == duration
                and r["kind"] == "r2"
                and r["repeat"] == repeat
                and r.get("summary")
            ),
            None,
        )
        if not protective:
            protective = next(
                (
                    r
                    for r in old
                    if r["case"] == case
                    and r["duration"] == duration
                    and r["kind"] == "nominal"
                    and r.get("repeat") == 0
                ),
                None,
            )
        rows.append(
            dict(
                case=case,
                duration_s=duration,
                method=kind,
                repeat=repeat,
                run=entry.get("run"),
                run_status=current,
                continuous_motion=dict(
                    status=numerical_status,
                    reason="Frozen accepted-state KKT/resistance, pair/WALL sweeps and zero projection; full requested duration required"
                    if kind != "reference"
                    else "Legacy benchmark has different projection/contact audit semantics",
                    state_integrity=state_audit.get("state_integrity") if state_audit else None,
                ),
                all_active_accuracy=dict(
                    status="NOT_EVALUABLE",
                    reason="Independent global trajectory error bound not established",
                    comparisons=[r["samples"] for r in comp],
                    comparison_scope="Only comparisons containing this exact run; an empty list means no completed comparison for this run",
                ),
                cluster_accuracy=dict(
                    status="NOT_EVALUABLE",
                    reason="Reference consistency is recorded separately from truth",
                    reference_consistency=[
                        dict(name=r["name"], duration_s=r["duration_s"], groups=r["cluster_consistency"])
                        for r in comp
                    ],
                ),
                reference_consistency=dict(
                    status="FAIL",
                    label="BASELINE_DISAGREEMENT",
                    reason="Archived/fresh lineage discrepancy remains; no reference selected to favor R3",
                ),
                performance_vs_frozen=ratio(ref, entry, duration)
                if kind == "nominal"
                else dict(
                    status="NOT_EVALUABLE", speedup=None, reason="This row is not a nominal timing candidate"
                ),
                performance_vs_R2=ratio(protective, entry, duration)
                if kind == "nominal"
                else dict(
                    status="NOT_EVALUABLE", speedup=None, reason="This row is not a nominal timing candidate"
                ),
                promotion=dict(
                    status="FAIL" if kind in ["nominal", "refined"] else "NOT_RUN",
                    reason="Independent all-active accuracy unavailable; numerical completion alone does not authorize promotion",
                ),
                accepted_summary=summary,
                charged_wall_s=entry.get("charged_wall_s"),
            )
        )
    performance = []
    for case in ["B2", "B3"]:
        new = [
            r
            for r in ledger
            if r["case"] == case
            and r["kind"] == "nominal"
            and r["duration_s"] == 0.005
            and completed(r, 0.005)
        ]
        ref = [
            r
            for r in ledger
            if r["case"] == case
            and r["kind"] == "reference"
            and r["duration_s"] == 0.005
            and completed(r, 0.005)
        ]
        ratios = [ratio(next((x for x in ref if x["repeat"] == r["repeat"]), None), r, 0.005) for r in new]
        full = (
            len(new) >= 3
            and len(ref) >= 3
            and len(ratios) >= 3
            and all(r["speedup"] is not None for r in ratios)
        )
        performance.append(
            dict(
                case=case,
                duration_s=0.005,
                complete_new_repeats=len(new),
                complete_reference_repeats=len(ref),
                new_wall_s=[r["charged_wall_s"] for r in new],
                reference_wall_s=[r["charged_wall_s"] for r in ref],
                paired=ratios,
                status="PASS"
                if full and all(r["speedup"] >= 5 for r in ratios)
                else "FAIL"
                if full
                else "INCOMPLETE",
                median_completion_speedup=float(
                    np.median([r["charged_wall_s"] for r in ref])
                    / np.median([r["charged_wall_s"] for r in new])
                )
                if full
                else None,
                scope="At least three completed nominal and complete-reference repeats required; incomplete references excluded",
            )
        )
    recovered = all(
        any(
            r["case"] == case
            and r["method"] == "nominal"
            and r["duration_s"] == 0.005
            and r["continuous_motion"]["status"] == "PASS"
            for r in rows
        )
        for case in ["B2", "B3"]
    )
    outcome = []
    if recovered:
        outcome.append("CONTINUATION_RECOVERED_ACCURACY_UNRESOLVED")
        if any(r["status"] != "PASS" for r in performance):
            outcome.append("NUMERICAL_RECOVERED_PERFORMANCE_PENDING")
    diagnostics = []
    for window in ["initial", "stable", "S0", "event", "B3_global_costly", "B3_serial_costly"]:
        path = BASE / "diagnosis" / f"DYNAMIC_INSTRUMENTED_{window}.json"
        d = read(path) if path.exists() else None
        diagnostics.append(
            dict(
                window=window,
                status="NOT_RUN"
                if d is None
                else "PASS"
                if d["summary"]["status"] == "COMPLETE"
                else "INCOMPLETE",
                source=str(path) if d else None,
                sha256=sha(path) if d else None,
                reason="Post-campaign instrumented evidence is separate from formal completion timing",
            )
        )
    pipeline_path = BASE / "diagnosis/finish_pipeline/RESULT.json"
    pipeline = read(pipeline_path) if pipeline_path.exists() else None
    pipeline_closed = bool(pipeline and pipeline["complete"])
    result = dict(
        campaign_complete=campaign["complete"],
        delivery_state="FINAL" if campaign["complete"] and pipeline_closed else "INTERIM",
        post_campaign_diagnostics=diagnostics,
        post_campaign_pipeline_closed=pipeline_closed,
        outcome=outcome or ["INCOMPLETE"],
        promotion=dict(status="FAIL", reason="Independent all-active global accuracy unavailable"),
        rows=rows,
        performance=performance,
        protected_R2_repeat_performance=protected_performance(ledger),
        tests=dict(status=tests["status"], counts=tests["summary"], source_bound=True),
        holdout_status=holdout["status"],
        holdout_scope="Six preselected full-system numerical states; trajectory-level lineage was previously inspected",
        spent_formal_s=campaign["spent_s"],
        protocol_sha256=sha(protocol_path()),
        preserved_failures="R2 nominal/refined B2 failures; R3 finest initial local budgets; strict bitwise restart FAIL with separate 1-ULP resolution audit; early loaded-source correction; original failed test invocation",
        no_extrapolated_speedup=True,
        no_automatic_promotion=True,
    )
    write(BASE / "final/MOTION_SOLVER_R3_GATE.json", result)
    print(dict(outcome=result["outcome"], complete=campaign["complete"], performance=performance), flush=True)


def protected_performance(ledger):
    rows = []
    for case in ["B2", "B3"]:
        new = [
            r
            for r in ledger
            if r["case"] == case
            and r["kind"] == "nominal"
            and r["duration_s"] == 0.005
            and completed(r, 0.005)
        ]
        old = [
            r
            for r in ledger
            if r["case"] == case and r["kind"] == "r2" and r["duration_s"] == 0.005 and completed(r, 0.005)
        ]
        paired = [ratio(next((o for o in old if o["repeat"] == n["repeat"]), None), n, 0.005) for n in new]
        complete = len(new) >= 3 and len(old) >= 3 and all(r["speedup"] is not None for r in paired)
        rows.append(
            dict(
                case=case,
                duration_s=0.005,
                new_wall_s=[n["charged_wall_s"] for n in new],
                protected_R2_wall_s=[o["charged_wall_s"] for o in old],
                paired=paired,
                complete_new_repeats=len(new),
                complete_R2_repeats=len(old),
                status="PASS"
                if complete and all(r["speedup"] >= 5 for r in paired)
                else "FAIL"
                if complete
                else "NOT_EVALUABLE",
                median_completion_speedup=float(
                    np.median([o["charged_wall_s"] for o in old])
                    / np.median([n["charged_wall_s"] for n in new])
                )
                if complete
                else None,
                reason="Three complete equal-start equal-duration R2 and R3 repeats required; historical incomplete B2 R2 is not extrapolated",
            )
        )
    return rows


if __name__ == "__main__":
    main()

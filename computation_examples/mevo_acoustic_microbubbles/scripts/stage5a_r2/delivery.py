"""Post-run evidence aggregation. Never mutates solver results or acceptance gates."""

from collections import Counter
from datetime import datetime, timezone
import html
import json
import numpy as np
import h5py
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage5a_r2.validation_figures import local
from scripts.stage4l.dense_probe import write, sha


def read(path, default=None):
    path = local(path)
    return json.loads(path.read_text()) if path.exists() else default


def lines(path):
    path = local(path)
    return [json.loads(s) for s in path.read_text().splitlines() if s] if path.exists() else []


def stats(values):
    a = np.asarray([v for v in values if v is not None], dtype=float)
    assert np.isfinite(a).all()
    return (
        dict(
            n=len(a),
            minimum=float(a.min()),
            median=float(np.median(a)),
            maximum=float(a.max()),
            mean=float(a.mean()),
            std_sample=float(a.std(ddof=1)) if len(a) > 1 else None,
        )
        if len(a)
        else None
    )


def numerical(s):
    return bool(
        s
        and s["status"] == "COMPLETE"
        and s["accepted_steps"] > 0
        and s["accepted_kkt_max"] <= 1e-9
        and s["accepted_linear_residual_max"] <= 1e-11
        and (s["accepted_nonlinear_kkt_max"] is None or s["accepted_nonlinear_kkt_max"] <= 1e-9)
        and s["accepted_pair_swept_gap_min_m"] >= -1e-11
        and s["accepted_wall_swept_gap_bound_min_m"] >= -1e-11
        and s["projection_max_m"] == 0
        and s["source_inputs_unchanged"]
    )


def run_audit(e):
    p = local(e["run"])
    s = e["summary"]
    a = lines(p / "ACCEPTED.jsonl")
    r = lines(p / "REJECTED.jsonl")
    if e["kind"] == "reference":
        return dict(
            run=e["run"],
            kind=e["kind"],
            case=e["case"],
            duration_s=e["duration"],
            summary=s,
            total_wall_s=e["wall_s_with_process_start"],
            reference_scope="Frozen Stage4L, same environment; not experimental truth.",
        )
    completed = read(p / "COMPLETED_TRIALS.json", [])
    resources = read(p / "RESOURCE_SAMPLES.json", [])
    identity = read(p / "IDENTITY.json")
    with h5py.File(p / "trajectories.h5") as h:
        t = h["time_s"][:]
        mask = h["active"][:]
        ids = h["particle_id"][:]
        rad = h["radii_m"][:]
        active_finite = all(np.isfinite(h["positions_m"][i][m]).all() for i, m in enumerate(mask))
        mask_constant = bool(np.all(mask == mask[0]))
    ledger_clock = bool(
        len(t) == len(a) + 1
        and np.all(np.diff(t) > 0)
        and np.allclose(t[1:], [q["time_s"] for q in a], rtol=0, atol=0)
        and abs(sum(q["dt_s"] for q in a) - (t[-1] - t[0])) < 2e-12
    )
    accepted = dict(
        count=len(a),
        numerical_qp_kkt=stats(q["kkt"] for q in a),
        original_S1_endpoint_kkt=stats(q["nonlinear_kkt"] for q in a),
        original_S0_start_kkt=stats(q["kkt"] for q in a if q["substep_method"] == "S0"),
        force_residual_l2_N=stats(q["force_residual_l2_N"] for q in a),
        resistance_relative=stats(q["linear_residual"] for q in a),
        dt_s=stats(q["dt_s"] for q in a),
        method_counts=dict(Counter(q["substep_method"] for q in a)),
    )
    rejected_completed = [q for q in r if q["reason"] == "LOCAL_ERROR"]
    rejected_stats = dict(
        count=len(rejected_completed),
        original_equilibrium_kkt=stats(
            q["diagnostics"].get("original_equilibrium_kkt") for q in rejected_completed
        ),
        numerical_qp_kkt=stats(q["diagnostics"].get("kkt") for q in rejected_completed),
        force_residual_l2_N=stats(
            q["diagnostics"].get("physical_force_residual_l2_N") for q in rejected_completed
        ),
        note="Only fully computed but rejected LOCAL_ERROR trials. Failed nonlinear trials are separate.",
    )
    failures = [
        dict(
            reason=q["reason"],
            time_s=q["time_s"],
            dt_s=q["dt_s"],
            trial_wall_s=q["trial_wall_s"],
            diagnostics=q["diagnostics"],
        )
        for q in r
        if q["reason"] != "LOCAL_ERROR"
    ]
    trial_cost = dict(
        completed_including_accepted_macro_and_rejected_local_error=len(completed) + len(rejected_completed),
        qp_iterations=sum(q.get("qp_iterations", 0) for q in completed)
        + sum(q["diagnostics"].get("qp_iterations", 0) for q in rejected_completed),
        nonlinear_iterations=sum(q.get("nonlinear_iterations", 0) for q in completed)
        + sum(q["diagnostics"].get("nonlinear_iterations", 0) for q in rejected_completed),
        accepted_substep_qp_iterations=sum(q["qp_iterations"] for q in a),
        accepted_substep_native_fallbacks=sum(q["fallbacks"] for q in a),
        accepted_substep_workspace_hits=sum(q["workspace_hits"] for q in a),
        accepted_substep_component_solves=sum(q["component_solves"] for q in a),
        accepted_substep_workspace_hit_fraction=sum(q["workspace_hits"] for q in a)
        / max(1, sum(q["component_solves"] for q in a)),
        accepted_substep_native_timing_s={
            key: sum(q["qp_timings"][key] for q in a)
            for key in ["setup_s", "update_s", "solve_s", "polish_s"]
        },
        accepted_component_size_statistics=stats(size for q in a for size in q["component_sizes"]),
        total_constrained_stage_calls=s["qp_calls"],
        limitation="Iteration totals cover completed trials. Native timing/workspace/component subcounts cover adopted substeps only, not all discarded full-step trials. Total qp timer and wall time include every attempt. Terminal/nonconverged trials keep available diagnostics; unavailable internal factorization/iteration histories are not invented. polish_s is external polish; native OSQP polish is included in solve_s.",
    )
    return dict(
        run=e["run"],
        case=e["case"],
        kind=e["kind"],
        duration_s=e["duration"],
        summary=s,
        total_wall_s=e["wall_s_with_process_start"],
        process_start_and_outer_wrapper_s=e["wall_s_with_process_start"] - s["wall_s"],
        identity=identity,
        accepted_state=accepted,
        completed_rejected_trials=rejected_stats,
        failed_trial_candidates=failures,
        failed_trial_candidate_statistics=dict(
            count=len(failures),
            last_original_equation_kkt=stats(
                q["diagnostics"].get("terminal_original_nonlinear_kkt") for q in failures
            ),
            last_physical_force_residual_l2_N=stats(
                q["diagnostics"].get("terminal_physical_force_residual_l2_N") for q in failures
            ),
        ),
        terminal_run_failure=s["failure"],
        terminal_run_error=s.get("error"),
        terminal_status=s["status"],
        last_rejected_reason=r[-1]["reason"] if r else None,
        failed_candidate_residual_limitations="Sweep-rejected exceptions preserve their continuous geometry bounds and dt, but do not retain native-candidate KKT/force fields. Those residuals remain unavailable, not zero. Nonlinear failures retain their recorded final original-equation residuals.",
        trial_cost=trial_cost,
        geometry=dict(
            finite_active_states=bool(active_finite),
            constant_IDs=ids.tolist(),
            constant_radii_m=rad.tolist(),
            active_count=int(mask[0].sum()),
            active_mask_constant=mask_constant,
            physical_time_ledger_consistent=ledger_clock,
            continuous_pair_and_wall_certificates=all(q["safety"]["pass_gate"] for q in a) if a else None,
        ),
        resource_statistics={
            k: stats(q.get(k) for q in resources)
            for k in [
                "process_cpu_percent",
                "host_cpu_percent",
                "gpu_percent",
                "rss_bytes",
                "gpu_memory_MiB",
                "active_workers",
            ]
        },
        resource_sampling_note="Raw first CPU sample follows immediate psutil priming (very short interval); an initial spike is not sustained parallel usage. Report median/whole-process CPU time alongside extrema. Host/GPU utilization is device-wide.",
        note="Inactive NaN sentinel positions are excluded from finite-active-state check; IDs/radii/active masks are retained.",
    )


def performance(entries):
    rows = []
    for d in sorted({e["duration"] for e in entries}):
        for c in ["B2", "B3"]:
            times = {
                k: [
                    e["wall_s_with_process_start"]
                    for e in entries
                    if e["case"] == c
                    and e["duration"] == d
                    and e["kind"] == k
                    and e.get("summary", {}).get("status") == "COMPLETE"
                ]
                for k in ["nominal", "reference"]
            }
            comparable = all(times.values())
            speed = float(np.median(times["reference"]) / np.median(times["nominal"])) if comparable else None
            rows.append(
                dict(
                    case=c,
                    duration_s=d,
                    nominal_wall_s=times["nominal"],
                    reference_wall_s=times["reference"],
                    nominal_stats=stats(times["nominal"]),
                    reference_stats=stats(times["reference"]),
                    completion_speedup=speed,
                    repeatability_sufficient=all(len(v) >= 3 for v in times.values()),
                    status="PASS" if comparable and speed >= 5 else "FAIL" if comparable else "INCOMPLETE",
                    scope="Same fresh environment, BLAS=1, CuPy pair backend, serial jobs; all starts, retries and IO included. Historical times excluded.",
                )
            )
    return rows


def main():
    campaign_path = sorted((BASE / "final").glob("CAMPAIGN_V*/RESULT.json"))[-1]
    version = campaign_path.parent.name.rsplit("_", 1)[-1]
    campaign = read(campaign_path)
    assert campaign and campaign["complete"]
    entries = [e for e in campaign["ledger"] if e.get("summary")]
    audits = [run_audit(e) for e in entries]
    write(BASE / "final/RUN_AUDIT.json", audits)
    perf = performance(entries)
    write(BASE / "final/PERFORMANCE.json", dict(rows=perf, budget_spent_s=campaign["spent_s"]))
    partial = read(BASE / "final/PARTIAL_COMPARISONS.json")
    repetitions = read(BASE / "final/REPEAT_AUDIT.json")
    prefix_disagreement = read(BASE / "final/PARTIAL_BASELINE_DISAGREEMENT.json")
    stages = {}
    blockers = []
    endpoint = read(BASE / "final/INDEPENDENT_ENDPOINT_AUDIT.json")
    endpoint_by_run = {r["run"]: r for r in endpoint["rows"]} if endpoint else {}
    for duration in [0.001, 0.005]:
        raw = campaign["stages"].get(str(duration))
        if not raw:
            stages[str(duration)] = dict(
                continuous_motion_status="NOT_RUN", scientific_status="NOT_EVALUABLE", cases={}
            )
            blockers.append(f"{duration * 1e3:g} ms stage not run: prerequisite gates not passed")
            continue
        cases = {}
        for case in ["B2", "B3"]:
            check = raw["cases"][case]
            num = check["nominal_numerical"] and check["refined_numerical"]
            members = [
                r
                for r in audits
                if r["case"] == case and r["duration_s"] == duration and r["kind"] != "reference"
            ]
            evidence_ok = all(
                r["geometry"]["finite_active_states"]
                and r["geometry"]["active_mask_constant"]
                and r["geometry"]["physical_time_ledger_consistent"]
                and r["geometry"]["continuous_pair_and_wall_certificates"]
                for r in members
            )
            failed = any(r["summary"]["status"] == "FAIL" for r in members) or not evidence_ok
            failed |= any(
                not endpoint_by_run[r["run"]]["pass_gate"] for r in members if r["run"] in endpoint_by_run
            )
            motion_status = "FAIL" if failed else "PASS" if num else "INCOMPLETE"
            science_keys = ["self_convergence", "legacy_comparison", "refined_legacy_comparison"]
            reports = [check.get(k) for k in science_keys]
            sci = (
                "PASS"
                if all(r and r["pass_scientific"] for r in reports) and check.get("actual_grid_refined")
                else (
                    "FAIL"
                    if any(r and r["status"] == "FAIL" for r in reports)
                    else "INCOMPLETE"
                    if any(reports)
                    else "NOT_EVALUABLE"
                )
            )
            repeat_rows = (
                [
                    r
                    for r in repetitions["rows"]
                    if r["case"] == case and abs(r["compared_duration_s"] - duration) < 1e-12
                ]
                if repetitions
                else []
            )
            repeated_science = (
                (
                    "PASS"
                    if all(r["pass_scientific"] for r in repeat_rows)
                    else "FAIL"
                    if any(r["status"] == "FAIL" for r in repeat_rows)
                    else "INCOMPLETE"
                )
                if repeat_rows
                else "NOT_RUN"
            )
            if repeated_science == "FAIL":
                sci = "FAIL"
            elif repeated_science == "INCOMPLETE" and sci == "PASS":
                sci = "INCOMPLETE"
            cases[case] = dict(
                continuous_motion_status=motion_status,
                scientific_status=sci,
                interpretation=check["status"],
                interpretation_scope="Original frozen campaign status; additional post-timing audits may add a delivery failure without rewriting that campaign.",
                actual_grid_refined=check.get("actual_grid_refined"),
                comparisons={k: check.get(k) for k in science_keys},
                repeated_scientific_status=repeated_science,
            )
            if motion_status != "PASS" or sci != "PASS":
                blockers.append(
                    f"{case} {duration * 1e3:g} ms: numerical={num}, science={sci}, frozen campaign={check['status']}, repeated science={repeated_science}"
                )
        stages[str(duration)] = dict(
            cases=cases,
            continuous_motion_status="PASS"
            if all(c["continuous_motion_status"] == "PASS" for c in cases.values())
            else "FAIL"
            if any(c["continuous_motion_status"] == "FAIL" for c in cases.values())
            else "INCOMPLETE",
            scientific_status="PASS"
            if all(c["scientific_status"] == "PASS" for c in cases.values())
            else "FAIL"
            if any(c["scientific_status"] == "FAIL" for c in cases.values())
            else "INCOMPLETE",
        )
    tests = read(BASE / f"final/TEST_GATE_{version}.json")
    holdout = read(BASE / f"final/HOLDOUT_{version}.json")
    manufactured = read(BASE / "diagnosis/manufactured_final/GATE.json")
    five = stages["0.005"]
    p5 = [r for r in perf if r["duration_s"] == 0.005]
    pstatus = (
        "PASS"
        if len(p5) == 2 and all(r["status"] == "PASS" for r in p5)
        else ("FAIL" if any(r["status"] == "FAIL" for r in p5) else "INCOMPLETE")
    )
    if pstatus != "PASS":
        blockers.append("Both 5 ms comparable >=5x timing gates are not established")
    if not endpoint or endpoint["status"] != "PASS":
        blockers.append("Independent endpoint audit incomplete or failed")
    if not manufactured or manufactured["status"] != "PASS":
        blockers.append("Final-kernel manufactured convergence audit incomplete or failed")
    if not repetitions or repetitions["status"] != "PASS":
        blockers.append("Post-timing scientific repeatability audit incomplete or failed")
    if prefix_disagreement:
        for row in prefix_disagreement["rows"]:
            if row["interpretation"] == "BASELINE_DISAGREEMENT":
                blockers.append(
                    f"{row['case']} BASELINE_DISAGREEMENT on the actual {row['compared_duration_s'] * 1e3:.8g} ms common prefix despite same-clock refinement; separate from uncompleted 5 ms endpoint gates"
                )
    motion_pass = five["continuous_motion_status"] == "PASS"
    science_pass = five["scientific_status"] == "PASS"
    promoted = (
        motion_pass
        and science_pass
        and pstatus == "PASS"
        and endpoint
        and endpoint["status"] == "PASS"
        and holdout["status"] == "PASS"
        and tests.get("status") == "PASS"
        and manufactured
        and manufactured["status"] == "PASS"
        and repetitions
        and repetitions["status"] == "PASS"
    )
    native = read(BASE / "diagnosis/native_v1/SUMMARY.json")
    fixed = sorted(
        {
            r["fixture"]
            for r in native
            if r["solver"] == "prox_dense" and r["eps"] == 1e-12 and r["kkt"] <= 1e-9
        }
    )
    promoted = bool(promoted and len(fixed) == 3)
    identity = read(BASE / f"final/VALIDATION_PROTOCOL_{version}.json")
    gate = dict(
        schema_version=1,
        time_utc=datetime.now(timezone.utc).isoformat(),
        implementation_status="PASS" if tests.get("status") == "PASS" else "FAIL",
        qp_regression_status="PASS" if len(fixed) == 3 else "FAIL",
        continuous_motion_status=five["continuous_motion_status"],
        scientific_status=five["scientific_status"],
        performance_status=pstatus,
        promotion="PASS" if promoted else "FAIL",
        promoted=bool(promoted),
        outcome="PROMOTED"
        if promoted
        else "NUMERICAL_PASS_PERFORMANCE_PENDING"
        if motion_pass and science_pass
        else "VALIDATION_INCOMPLETE_OR_FAILED",
        blocking_reasons=blockers,
        stages=stages,
        tests=tests,
        holdout_status=holdout["status"],
        independent_endpoint_audit_status=endpoint["status"] if endpoint else "NOT_RUN",
        manufactured_convergence_status=manufactured["status"] if manufactured else "NOT_RUN",
        repeated_scientific_status=repetitions["status"] if repetitions else "NOT_RUN",
        performance=perf,
        partial_prefix_comparisons=partial,
        partial_prefix_baseline_disagreement=prefix_disagreement,
        archived_full_duration_science=read(BASE / "final/ARCHIVED_SCIENCE.json"),
        timing_limitations=[
            f"{r['case']} {r['duration_s'] * 1e3:g} ms has fewer than three complete repeats per method"
            for r in perf
            if not r["repeatability_sufficient"]
        ],
        executed_source_commit=(ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip(),
        executed_source_commits=sorted({r["identity"]["source_commit"] for r in audits if "identity" in r}),
        campaign=str(campaign_path),
        validation_version=version,
        final_archive_identity="Resolve the Git commit containing this artifact; execution identity is separate. See DELIVERY_IDENTITY.json for archive parent.",
        frozen_protocol_sha256=sha(BASE / f"final/VALIDATION_PROTOCOL_{version}.json"),
        protocol=identity,
        scope="Stages refer to original B2/B3 t=3 s states, active 186/189; hybrid S1/S0 is explicit in each ledger. Computational gates do not validate dense secondary physics.",
        policy="PASS, FAIL, NOT_RUN, NOT_EVALUABLE and INCOMPLETE retain separate meanings; no projection, no tolerance relaxation, no extrapolated speedup.",
    )
    write(BASE / "final/MOTION_SOLVER_R2_GATE.json", gate)
    report(gate, entries)
    index()
    print(
        json.dumps(
            {
                k: gate[k]
                for k in [
                    "implementation_status",
                    "qp_regression_status",
                    "continuous_motion_status",
                    "scientific_status",
                    "performance_status",
                    "promotion",
                    "outcome",
                ]
            },
            indent=2,
        )
    )


def report(gate, entries):
    text = [
        "# Stage 5A-R2 measured outcome",
        "",
        f"Outcome: **{gate['outcome']}**. Promotion: **{gate['promotion']}**.",
        "",
        f"Executed source: `{gate['executed_source_commit']}`. Final archive commit is the commit containing this report; it is not the executed source identity.",
        "",
        f"Fresh bound tests: {gate['tests']['passed']} passed, {gate['tests']['skipped_or_xfail']} retained old xfail. Original-QP regressions: {gate['qp_regression_status']}; holdout: {gate['holdout_status']}; independent final-equation audit: {gate['independent_endpoint_audit_status']}; final-kernel manufactured convergence: {gate['manufactured_convergence_status']}.",
        "",
        "Native packages: proxsuite/ProxQP 0.7.2, OSQP 1.1.3, Clarabel 0.11.1; float64. The selected solver is bounded component ProxQP, with certified native alternatives. Internal status never replaces the fixed original KKT check.",
        "",
        "The final method uses explicit acoustic/flow forcing and endpoint implicit resistance (S1), with bounded start-state S0 recovery at nonsmooth near-wall normal transitions. Every adopted substep has step-doubling and continuous pair/WALL checks. The old tolerance remains 1e-11 m; no projection is applied.",
        "",
        "| Case | Kind | Target (ms) | Completed (ms) | Accepted / rejected | Min / median dt (s) | Total wall (s) | Status |",
        "|---|---|---:|---:|---:|---|---:|---|",
    ]
    for e in entries:
        s = e["summary"]
        dt = s.get("actual_dt_quantiles_s")
        dts = f"{dt[0]:.5g} / {dt[1]:.5g}" if dt else "See raw ledger"
        text.append(
            f"| {e['case']} | {e['kind']} repeat {e['repeat']} | {e['duration'] * 1e3:g} | {s.get('completed_duration_s', s.get('elapsed_s', 0)) * 1e3:.8g} | {s.get('accepted_steps', '—')} / {s.get('rejected_trials', '—')} | {dts} | {e['wall_s_with_process_start']:.3f} | {s['status']} |"
        )
    text += [
        "",
        "## Accepted-state numerical evidence",
        "",
        "| Case / target / kind | QP / original S1 KKT max | Resistance residual max | Force residual max (N) | Pair / WALL swept minimum (m) | Projection (m) |",
        "|---|---|---|---|---|---|",
    ]

    def fmt(v):
        return f"{v:.6g}" if v is not None else "null / not applicable"

    for e in entries:
        if e["kind"] == "reference" or e["repeat"] != 0:
            continue
        s = e["summary"]
        text.append(
            f"| {e['case']} / {e['duration'] * 1e3:g} ms / {e['kind']} | {fmt(s['accepted_kkt_max'])} / {fmt(s['accepted_nonlinear_kkt_max'])} | {fmt(s['accepted_linear_residual_max'])} | {fmt(s['accepted_force_residual_l2_max_N'])} | {fmt(s['accepted_pair_swept_gap_min_m'])} / {fmt(s['accepted_wall_swept_gap_bound_min_m'])} | {fmt(s['projection_max_m'])} |"
        )
    text += [
        "",
        "These are accepted states only. Complete LOCAL_ERROR rejections and terminal/nonconverged candidates are separately tabulated in RUN_AUDIT.json; failed candidates never contribute to accepted-state maxima. S0, if used, has a start-equation certificate and no endpoint nonlinear KKT. Raw negative gaps inside the unchanged 1e-11 m tolerance are retained.",
    ]
    text += [
        "",
        "## Scientific comparison",
        "",
        "| Duration / case | Comparison | Rg endpoint / series (%) | Q1 / Q2 / Q3 (%) | Centroid / budget (nm) | Status |",
        "|---|---|---|---|---|---|",
    ]
    for duration, stage in gate["stages"].items():
        for c, v in stage["cases"].items():
            for name, r in v["comparisons"].items():
                if not r:
                    continue
                text.append(
                    f"| {float(duration) * 1e3:g} ms / {c} | {name} | {r['Rg_endpoint_relative'] * 100:.5g} / {r['Rg_series_relative_L2'] * 100:.5g} | {' / '.join(f'{x * 100:.5g}' for x in r['pair_quartiles_relative'])} | {r['centroid_difference_m'] * 1e9:.5g} / {r['centroid_budget_m'] * 1e9:.5g} | {r['status']} |"
                )
    text += [
        "",
        "The original group membership and all particles are preserved. Rg and centroid have different units/thresholds. Comparisons use actual states at the common clock, with explicit sample-age bounds; no coordinate interpolation. Frozen Stage4L is a numerical reference, not experimental truth.",
        "Passing aggregate group metrics does not mean every individual trajectory is identical. The reports also retain maximum same-ID position difference, rigid-fit residual and neighbor-graph changes; R6/R7 display them without geometric smoothing.",
        "",
        "V3 ordinary B2 remains BASELINE_DISAGREEMENT (1 ms Q1 difference 2.20205%). V4/V4A use a documented tighter time-error configuration of the unchanged kernel. V4's rejected entry launches remain recorded; V4A fixes entry routing, with fresh bound tests and explicitly reused V4 holdouts. Earlier results are never relabelled as later configuration runs.",
        "Optional [partial-prefix comparisons](PARTIAL_COMPARISONS.json) describe actual common intervals when a 5 ms run stops early. They retain exact source frames and sampling uncertainty, and never upgrade a complete-duration gate.",
        "[Same-clock prefix refinement analysis](PARTIAL_BASELINE_DISAGREEMENT.json) checks whether a legacy disagreement persists while the new method converges on exactly that comparison clock, including conservative sampling bounds. BASELINE_DISAGREEMENT remains separate from uncompleted full-duration gates; neither reference truth nor new-method invalidity is inferred automatically.",
        "When a fresh reference cannot finish within budget, [archived full-duration science](ARCHIVED_SCIENCE.json) separately compares complete new 5 ms trajectories against the immutable original Stage4L data. These supplementary comparisons do not change the frozen campaign or use historical runtime as a speedup denominator.",
    ]
    for field, title in [
        ("partial_prefix_comparisons", "Supplementary actual common prefixes (not complete 5 ms gates)"),
        ("archived_full_duration_science", "Supplementary archived-reference science (not fresh timing)"),
    ]:
        extra = gate.get(field)
        if not extra or not extra["rows"]:
            continue
        text += [
            "",
            f"### {title}",
            "",
            "| Comparison | Actual common duration (ms) | Rg endpoint / series (%) | Q1 / Q2 / Q3 (%) | Centroid / budget (nm) | Scoped status |",
            "|---|---:|---|---|---|---|",
        ]
        for r in extra["rows"]:
            text.append(
                f"| {r['label']} | {r['compared_duration_s'] * 1e3:.8g} | {r['Rg_endpoint_relative'] * 100:.5g} / {r['Rg_series_relative_L2'] * 100:.5g} | {' / '.join(f'{x * 100:.5g}' for x in r['pair_quartiles_relative'])} | {r['centroid_difference_m'] * 1e9:.5g} / {r['centroid_budget_m'] * 1e9:.5g} | {r['status']} |"
            )
    text += [
        "",
        "## Performance",
        "",
        f"[Post-timing scientific repeatability](REPEAT_AUDIT.json): **{gate['repeated_scientific_status']}**. All measured repeats remain in the timing ledger. Adaptive step counts alone are not evidence of identical trajectories; repeated runs are compared against the first same-method run and, for new trajectories, the fresh first legacy reference under the unchanged scientific thresholds.",
        "",
    ]
    for r in gate["performance"]:
        text.append(
            f"- {r['case']} {r['duration_s'] * 1e3:g} ms: completion speedup {r['completion_speedup']}; new/reference repetitions {len(r['nominal_wall_s'])}/{len(r['reference_wall_s'])}; three-repeat sufficiency {r['repeatability_sufficient']}."
        )
    text += [
        "",
        "Only completed equal-duration runs under the same environment enter speedup. Total wall includes process startup, failed trials, nonlinear work, I/O and monitoring. Raw error/wall subcounters overlap; R8 uses a mutually exclusive subset plus the measured remainder. Native solvers do not expose complete factorization/inner iteration histories through all APIs.",
        "",
        "## Blockers",
        "",
    ]
    text += [f"- {r}" for r in gate["blocking_reasons"]] or ["- None within this numerical scope."]
    text += [
        "",
        "## Evidence and reproduction",
        "",
        "- [Machine-readable gate](MOTION_SOLVER_R2_GATE.json), [all run audits](RUN_AUDIT.json), [performance](PERFORMANCE.json).",
        "- [Figure and data index](../figures/index.md). Actual accepted/rejected ledgers, source snapshots and checkpoint files retain each run ID.",
        f"- Numerical code and environment are sealed by VALIDATION_PROTOCOL_{gate['validation_version']}.json and TEST_GATE_{gate['validation_version']}.json; earlier failed holdouts and campaigns remain unchanged.",
        f"- Formal reproduction uses the exact commands in CAMPAIGN_{gate['validation_version']}/RESULT.json, in a separate output root with the recorded scientific inputs and dependency lock. Explicitly reused runs keep their earlier commands, run IDs and source commits. Do not rerun the freeze writer over an existing protocol or overwrite immutable run directories.",
        "- Postprocessing: `PYTHONPATH=src:. python -m scripts.stage5a_r2.validation_figures`, `PYTHONPATH=src:. python -m scripts.stage5a_r2.motion_visuals --video --duration 0.001`, `PYTHONPATH=src:. python -m scripts.stage5a_r2.delivery`.",
        "- [Reproduction details](../../../../scripts/stage5a_r2/README.md), including exact source identities, external input files and the immutable-output policy.",
        "- Large scientific input fields, full HDF5 trajectories, source snapshots and MP4 exports are indexed by actual path, bytes and SHA256 in ARTIFACT_MANIFEST.json; they are not all embedded in Git.",
        "",
        "Known root causes and bounded revisions are documented in `docs/STAGE5A_R2_ROOT_CAUSE_ANALYSIS.md` and `docs/STAGE5A_R2_MOTION_DEVELOPMENT.md`. Passing numerical checks would not validate the dense classical secondary-force model or a physical aggregation mechanism.",
        "",
    ]
    (BASE / "final/REPORT.md").write_text("\n".join(text))


def index():
    out = BASE / "figures"
    rows = []
    for p in sorted(out.glob("R*.png")):
        metadata = out / (p.stem + "_SOURCES.json")
        if not metadata.exists():
            metadata = out / (p.stem.split("_")[0][:2] + "_DATA_SOURCES.json")
        if not metadata.exists() and p.stem.startswith("R7"):
            metadata = out / "R7_DATA_SOURCES.json"
        assert metadata.exists(), f"Missing figure sources for {p.name}"
        script = (
            "diagnostic_figures.py"
            if p.stem.startswith(("R3", "R4"))
            else "motion_visuals.py"
            if p.stem.startswith("R7")
            else "validation_figures.py"
        )
        rows.append((p, metadata, script))
    md = [
        "# Stage 5A-R2 figures and actual motion",
        "",
        "[Final measured report](../final/REPORT.md) · [Separated gates](../final/MOTION_SOLVER_R2_GATE.json)",
        "",
        "All coordinates and statistics use actual accepted states. Raw sources and SHA256 are linked for each figure; source/run identities are not inferred from visual appearance.",
        "",
    ]
    body = [
        '<h1>Stage 5A-R2: evidence index</h1><p><a href="../final/REPORT.md">Measured report</a> · <a href="../final/MOTION_SOLVER_R2_GATE.json">Separated gates</a></p>'
    ]
    for p, metadata, script in rows:
        link = "../../../../scripts/stage5a_r2/" + script
        md.extend(
            [
                f"## {p.stem}",
                "",
                f"[PNG]({p.name}) · [PDF]({p.with_suffix('.pdf').name}) · [raw sources and definitions]({metadata.name}) · [generator]({link})",
                "",
                f"![{p.stem}]({p.name})",
                "",
            ]
        )
        body.append(
            f'<section><h2>{html.escape(p.stem)}</h2><p><a href="{p.with_suffix(".pdf").name}">PDF</a> · <a href="{metadata.name}">Raw sources and definitions</a> · <a href="{link}">Generator</a></p><img loading="lazy" src="{p.name}" alt="{p.stem}"></section>'
        )
    video = out / "R7_OLD_NEW_ACTUAL_MOTION_60fps.mp4"
    if video.exists():
        md += [
            "## Same-clock actual motion",
            "",
            f"[60 fps MP4]({video.name}) · [actual frame data](R7_VIDEO_ACTUAL_FRAMES.npz) · [sampling and camera provenance](R7_DATA_SOURCES.json)",
            "",
            "[Decoded-frame and coordinate audit](R7_VIDEO_VALIDATION.json). The MP4 and full frame archive are external exports when excluded by the size policy; local links work in the complete export folder. For GitHub readers, retrieve the actual paths, sizes and SHA256 from the [artifact manifest](../final/ARTIFACT_MANIFEST.json).",
            "",
        ]
        body.append(
            f'<h2>Same-clock actual motion</h2><video controls preload="metadata" src="{video.name}"></video>'
            "<p>The video is an external export when absent from the Git checkout. "
            '<a href="../final/ARTIFACT_MANIFEST.json">Actual file paths and SHA256</a> · '
            '<a href="R7_VIDEO_VALIDATION.json">Decoded-frame and coordinate audit</a></p>'
        )
    (out / "index.md").write_text("\n".join(md))
    (out / "index.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><title>Stage 5A-R2 evidence</title><style>body{max-width:1250px;margin:30px auto;font:16px system-ui;color:#26313c}img,video{width:100%}section{margin:3em 0}a{color:#275b89}</style><body>'
        + "".join(body)
        + "</body></html>"
    )


if __name__ == "__main__":
    main()

"""Generate the final report from closed measured ledgers, never predictions."""

import subprocess
from scripts.stage5a_r3.common import ROOT, BASE, read
from scripts.stage5a_r3.current import campaign


def value(x, scale=1, digits=6):
    return "unavailable" if x is None else f"{x * scale:.{digits}g}"


def main():
    c = campaign()
    gate = read(BASE / "final/MOTION_SOLVER_R3_GATE.json")
    assert c["complete"] and gate["campaign_complete"]
    assert gate["delivery_state"] == "FINAL", "Post-campaign diagnostics must be accounted for"
    comparisons = read(BASE / "final/FORMAL_COMPARISONS.json")["rows"]
    usage = read(BASE / "final/FORMAL_RUN_AUDIT.json")["rows"]
    resources = read(BASE / "final/RESOURCE_LAYOUT_BENCHMARK.json")
    budget = read(BASE / "final/BUDGET_LEDGER.json")
    current_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    lines = [
        "# Stage 5A-R3 final development and validation report",
        "",
        "Outcome: **" + " / ".join(gate["outcome"]) + "**. Promotion: **FAIL**.",
        "",
        "B2 continuous motion has been repaired without changing forces, radii, the WALL mesh, lubrication coefficients or tolerances. Independent global all-active accuracy remains NOT_EVALUABLE. The measured performance verdict is separate from numerical completion.",
        "",
        "## Identity and scope",
        "",
        "- Branch: `dev/stage5a-r3-wall-accuracy-performance`.",
        "- Actual R2 archive base: `07db8ed35f3a460a4ba315e42b15d257f8518b0a`.",
        "- V1 executable: `fa2fba5608037ed45143c8ede501cadf415e9180`.",
        "- V2 executable: `b5079ad19be81ffe9e3d4f5b7c548717b3d43d67`.",
        f"- Report generation source: `{current_sha}`; resolve the containing Git commit for the final archive SHA.",
        "- Server: `/workspace/mevo_acoustic_microbubbles_stage5a_r3`; immutable executable trees are `.execution/<SHA>`.",
        "- Results: `results/stage5/motion_solver_r3`; no old branch, checkpoint, failed gate or stopped service is resumed or overwritten.",
        "",
        "V2 changes B3 QP resource organization to component serial, retaining B2 global QP. BLAS and QP workers remain one, and CuPy supplies the pair force. Equivalent B2 V1 executions retain their original identity and cost; no budget reset or relabeling occurred. Formal runs begin at the original 3 s state and extend by at most 5 ms.",
        "",
        "## Actual completion and accepted numerical safety",
        "",
        "| Case | Duration (ms) | Method | Repeat | Status | Process-inclusive wall (s) | QP KKT max | Pair / WALL minimum (m) |",
        "|---|---:|---|---:|---|---:|---:|---|",
    ]
    for r in c["ledger"]:
        s = r.get("summary") or {}
        lines.append(
            f"| {r['case']} | {r['duration_s'] * 1000:g} | {r['kind']} | {r['repeat']} | {r['status']} | {value(r.get('charged_wall_s'))} | {value(s.get('accepted_kkt_max'))} | {value(s.get('accepted_pair_swept_gap_min_m'))} / {value(s.get('accepted_wall_swept_gap_bound_min_m'))} |"
        )
    lines += [
        "",
        "A missing process-inclusive time is not replaced by body timing. In particular, B2 V1 refined completed in 1152.959 s of measured run-body time while its parent scheduler was held. Body time plus 50 s is a conservative budget charge only. Legacy references use different projection/contact semantics, so absent R3 safety fields remain unavailable.",
        "",
        "S1 uses endpoint resistance with the frozen forcing/contact discretization. S0 uses start resistance and has no endpoint nonlinear KKT: NOT_APPLICABLE/null. QP KKT ≤1e−9, resistance residual ≤1e−11, pair/WALL tolerance 1e−11 m and zero projection remain unchanged. Rejected candidate residuals are separate from accepted maxima. HDF5 active IDs, finite coordinates, strict time increase and accepted-log correspondence are audited independently of the run summary.",
        "",
        "| Case / run | S0 accepted / duration (µs) | S1 accepted / duration (ms) | State integrity |",
        "|---|---:|---:|---|",
    ]
    for r in usage:
        methods = {m["method"]: m for m in r["methods"]}
        s0, s1 = methods["S0"], methods["S1"]
        lines.append(
            f"| `{r['run'].rsplit('/', 1)[-1]}` | {s0['accepted_substeps']} / {value(s0['physical_duration_s'], 1e6)} | {s1['accepted_substeps']} / {value(s1['physical_duration_s'], 1e3)} | {r['state_integrity']['status']} |"
        )
    lines += [
        "",
        "## Confirmed B2 failure and numerical repair",
        "",
        "The responsible bubble is ID 152. Its actual nearest WALL feature switches between finite faces 3178 and 13018, with a normal jump around 0.1343. Cold queries exactly match the stored-query gap/normal, refuting a stale-cache explanation. A 1.526 ps candidate is safe under the exact finite-triangle predicate but fails the old global roundoff guard; 1 ns and 50 ns candidates are genuinely unsafe. Therefore simply removing the guard is insufficient. The minimum-progress failure is an algorithmic 1 ps limit, not an unrepresentable timestamp.",
        "",
        "R3 rebuilds finite reachable hard-contact features at the actual position and uses directed-rounding separation bounds with exact binary64-rational segment/triangle fallback. Face/edge/vertex coverage, signs and open ports are tested. The nearest-wall lubrication model is unchanged; extra contact normals do not add extra wall resistance. No wall smoothing, projection, bubble deletion or synthetic accepted frames are used. See the root-cause report, G1–G3 and the full original failed ledgers in the review directory.",
        "",
        "## All-active precision and reference lineage",
        "",
        "Every comparison preserves original IDs and all 186/189 active particles, with the original 75/91-ID group and its complement also reported. Raw, translation-removed and rigid-fit errors are retained. Independent global trajectory bounds are NOT_ESTABLISHED, so these comparisons cannot promote a per-bubble reference.",
        "",
        "| Case / named comparison | Common duration (ms) | Endpoint raw RMS / max (nm) | Worst sampled ID / time (s) | Common-time bound available |",
        "|---|---:|---|---|---|",
    ]
    for r in comparisons:
        if not r.get("samples"):
            continue
        a = r["groups"]["all_active"]["errors_m"]["raw"]
        lines.append(
            f"| {r['case']} / {r['name']} | {r['duration_s'] * 1000:.6g} | {value(a['endpoint']['rms'], 1e9)} / {value(a['endpoint']['max'], 1e9)} | {a['worst_particle_id']} / {a['worst_time_s']:.9f} | {r.get('sampling_bound_certified', False)} |"
        )
    lines += [
        "",
        "B2 nominal/refined +5 ms endpoint RMS/max are 31.088/233.979 nm, exceeding frozen 17.8/89 nm targets. The first sampled maximum crossing is bracketed by +30..40 µs (ID 82); the worst sampled error is 669.608 nm at +70 µs (ID 211). This early event is distinct from the +4.31 ms WALL stall. Both selected IDs have no WALL face switch in +30..80 µs; no S0 or fallback occurs there.",
        "",
        "Six all-active local +20 µs replays use the actual nominal/refined +30 µs states at three precision levels. Adjacent finest RMS/max differences are 11.004/50.914 nm and 11.430/82.315 nm. Differences between the two actual starts grow from 4.888 nm RMS to 48.655 nm RMS under the finest local settings. Near-contact rearrangement is correlated with this sensitivity, but no unique global-error cause or preferred truth is established. Cold-cache/controller scope is explicit. See B2_EARLY_PRECISION_AUDIT.json and G5b.",
        "",
        "Archived Stage4L and fresh frozen Stage4L remain separate numerical lineages. Controlled tests identify early sensitivity to their 500/10 µs clock cutoff; same-input equation/QP checks agree to floating-point precision. This does not resolve the complete archived/fresh trajectory disagreement. At unmatched legacy times, a neighboring-frame chord is only a displacement scale unless a straight recorded path is established; consistency gates are withheld when that sampling bound is unavailable. Actual common stored frames are not interpolated.",
        "",
        "## Measured performance and resources",
        "",
        "| Case | Complete new / reference 5 ms repeats | New wall samples (s) | Reference wall samples (s) | Median completion speedup | Verdict |",
        "|---|---|---|---|---:|---|",
    ]
    for p in gate["performance"]:
        lines.append(
            f"| {p['case']} | {p['complete_new_repeats']} / {p['complete_reference_repeats']} | {', '.join(value(v) for v in p['new_wall_s'])} | {', '.join(value(v) for v in p['reference_wall_s'])} | {value(p['median_completion_speedup'])} | {p['status']} |"
        )
    lines += [
        "",
        "| Case | Complete protected R2 5 ms repeats | Protected R2 wall samples (s) | Median R2 / R3 speedup | Verdict |",
        "|---|---:|---|---:|---|",
    ]
    for p in gate["protected_R2_repeat_performance"]:
        lines.append(
            f"| {p['case']} | {p['complete_R2_repeats']} | {', '.join(value(v) for v in p['protected_R2_wall_s'])} | {value(p['median_completion_speedup'])} | {p['status']} |"
        )
    lines += [
        "",
        "Per-repeat comparisons against R2 and fresh Stage4L appear in MOTION_SOLVER_R3_GATE.json. Only complete equal-duration runs with verified identical initial arrays enter a speedup. Incomplete R2 B2 +5 ms is never extrapolated. Three repeats are required for the final performance gate; 5× is the target, not a promised result.",
        "",
        "The initial 50 ns layout benchmark did not represent B3's full-task cost. Complete 5 µs full/half1/half2 trials at +1.5/+2.5 ms gave global medians 1.618/1.321 s versus serial-component 0.351/0.288 s, with fixed-ID differences ≤2.43e−17 m. Persistent threaded component layouts were slower. The actual full tasks determine the final speedup; layout microbenchmarks do not.",
        "",
    ]
    for p in resources["pair_force"]:
        lines.append(
            f"- {p['case']} float64 GPU pair-force inclusive speedup: {p['force_only_speedup']:.4g}×. Transfers and synchronization are included; this is not trajectory speedup."
        )
    lines += [
        f"- Independent two-job throughput speedup: {resources['independent_throughput']['speedup']:.4g}×, separately labeled.",
        "- CPU affinity/cgroup, per-process CPU time, actual worker threads, memory and GPU transfer/kernel/synchronization records are retained. Device/host-wide utilization is not evidence of solver parallelism.",
        "- G7 and DYNAMIC_COST_PROFILE.json partition actual force/geometry/assembly/QP/nonlinear/full-half work without summing overlapping timers. Instrumented diagnostics are separate from formal completion timing.",
        "",
        "## Tests, budgets and provenance",
        "",
        f"Source-bound formal tests: **{gate['tests']['status']}**, counts `{gate['tests']['counts']}`. Six new numerical-state holdouts: **{gate['holdout_status']}**. Later reporting/sampling tests have a separate JUnit file; they do not relabel the frozen test executable.",
        "",
        f"Development conservative budget charge: {budget['development_charged_s']:.3f} s / 14400 s. Formal carried charges: `{budget['formal_charged_s']}`. Overall charged: {budget['total_charged_s']:.3f} s / 63000 s. Development allowances and the unmetered exploration reserve are budget accounting, not measured performance. No reset follows a failed candidate.",
        "",
        "The original strict bitwise restart failure (one coordinate ULP, 2.1684e−19 m) is retained alongside the separate resolution-based PASS and exact cache serialization PASS. Missing-fixture/import failed test invocations remain recorded. Two early local source ZIPs captured a later deployed driver after workers had loaded V1; original snapshots plus LOCAL_V1_PROVENANCE_CORRECTION.json and a reconstructed loaded-source archive remain. Later formal execution uses immutable trees.",
        "",
        "## Deliverables and remaining blockers",
        "",
        "- Required audit JSON and gate: `results/stage5/motion_solver_r3/final/`.",
        "- G1–G8 PNG/PDF, underlying actual data and source hashes: `figures/` and `figures/data/`.",
        "- Actual-radius event video and complete old/new videos, decoded timestamp/frame audits: `videos/`.",
        "- Full-system fixtures, candidate arrays, real WALL and complete failed-run ledgers: `review/`; standalone `Stage5A_R3_review.zip`.",
        "- Complete raw/server-only file paths and hashes: `final/ARTIFACT_MANIFEST.json`; exact execution/archive distinctions: `final/DELIVERY_IDENTITY.json`.",
        "- Reproduction commands: `scripts/stage5a_r3/README.md`; isolated dependency lock and license inventory retained.",
        "",
        "Remaining blockers are global all-active trajectory accuracy/reference disagreement and any failed measured performance gate. Numerical safety, global mass conservation, GPU agreement or a smooth animation does not establish local physical accuracy. No production/reference promotion and no longer-duration run is authorized by these partial successes.",
        "",
    ]
    (ROOT / "docs/STAGE5A_R3_FINAL_REPORT.md").write_text("\n".join(lines))
    print("Wrote docs/STAGE5A_R3_FINAL_REPORT.md from closed actual ledgers")


if __name__ == "__main__":
    main()

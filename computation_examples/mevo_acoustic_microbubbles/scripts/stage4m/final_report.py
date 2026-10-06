"""Finalize a gated Stage 4M outcome without inventing unrun formation data."""

from datetime import datetime, timezone
import json

import numpy as np

from scripts.stage4m.common import BASE, ROOT, write


def removal_summary(path):
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    t = np.array([r["time_s"] for r in rows])
    ratio = np.array(
        [np.nan if r["effective_classical_ratio"] is None else r["effective_classical_ratio"] for r in rows]
    )
    finite = np.isfinite(ratio)
    mean = float(np.sum(ratio[:-1] * np.diff(t)) / (t[-1] - t[0])) if len(t) > 1 and finite.all() else None
    return dict(
        saved_states=len(rows),
        initial=rows[0],
        final=rows[-1],
        ratio_min=float(np.min(ratio[finite])) if finite.any() else None,
        ratio_max=float(np.max(ratio[finite])) if finite.any() else None,
        time_weighted_left_hold_ratio=mean,
        all_pair_norms_finite=all(np.isfinite([r["S_classical_N"], r["S_effective_N"]]).all() for r in rows),
        semantics="Ratios use unordered pair-norm sums, not particle net force; left-hold time weighting is explicitly a diagnostic approximation.",
    )


def main():
    gate = json.loads((BASE / "comparisons/M0_GATE.json").read_text())
    if gate["continue_formation"]:
        raise RuntimeError("M0 passed: downstream formation/growth gates are required before a final report")
    data = json.loads((BASE / "comparisons/DENSE_COMPARISON.json").read_text())
    early = json.loads((BASE / "l0_vs_a0/LUBRICATION_ONLY_FORMATION_EFFECT.json").read_text())
    early_text = []
    for case, models in early["rows"].items():
        l0, a0 = models["L0"], models["A0"]
        early_text.append(
            f"{case}: largest cluster {l0['largest_cluster_size']} to {a0['largest_cluster_size']}; near-contact fraction {l0['near_contact_fraction']:.6f} to {a0['near_contact_fraction']:.6f}."
        )
    force, table = {}, []
    for case in ("B2", "B3"):
        force[case] = {}
        for policy in ("A0", "A1", "A2"):
            row = data[case][policy]
            s, ref = row["summary"], data[case]["A0"]["summary"]
            table.append(
                dict(
                    case=case,
                    policy=policy,
                    status=s["status"],
                    wall_s=s["wall_s"],
                    accepted_steps=s["accepted_steps"],
                    median_dt_s=s["median_dt_s"],
                    step_reduction=ref["accepted_steps"] / s["accepted_steps"]
                    if s["accepted_steps"]
                    else None,
                    wall_speedup=ref["wall_s"] / s["wall_s"],
                    max_rg_ratio=row.get("max_rg_ratio"),
                    contact_QP_calls=row.get("translation", {}).get("contact_solves"),
                    contact_QP_wall_share=row.get("translation", {}).get("total_contact_solve_time_s", 0)
                    / s["wall_s"],
                    projection_calls=s["projection"]["calls"],
                    multi_radius_projections=s["projection"]["multi_radius_corrections"],
                    max_projection_m=s["projection"]["max_correction_m"],
                )
            )
            if policy != "A0":
                force[case][policy] = removal_summary(
                    ROOT / s["run"] / "PAIR_FORCE_DIAGNOSTICS.jsonl"
                )
    write(BASE / "comparisons/SECONDARY_REMOVAL_SUMMARY.json", force)
    answers = dict(
        Q1=dict(
            status="EXISTING_EQUAL_TIME_COMPARISON_COMPLETE",
            answer="At 0.08 s, lubrication retains strong early aggregation. "
            + " ".join(early_text)
            + " Some formation metrics change, but a single checkpoint cannot establish a continuous formation-rate change.",
            evidence="l0_vs_a0/LUBRICATION_ONLY_FORMATION_EFFECT.json",
        ),
        Q2=dict(
            status="DENSE_DIAGNOSTIC_COMPLETE"
            if all(data[c]["A1"]["status"] == "COMPLETE" for c in data)
            else "DENSE_DIAGNOSTIC_INCOMPLETE",
            answer="A1 runtime and accepted-step changes are listed per case; it is a hard-switch diagnostic only, with no longer runs.",
            cases={c: next(r for r in table if r["case"] == c and r["policy"] == "A1") for c in data},
        ),
        Q3=dict(
            status="A2_RUNTIME_GATE_FAILED",
            answer="A2 did not meet the required fivefold continuation gate. Reduced pair-norm sums do not remove the observed small-step bottleneck sufficiently for a long campaign.",
            cases=gate["cases"],
        ),
        Q4=dict(
            status="NOT_RUN_M0_GATE_FAILED",
            answer="Whether A2 forms clusters from the original initial population remains unanswered. Legacy-formed dense checkpoints cannot answer formation.",
        ),
        Q5=dict(
            status="NOT_RUN_M0_GATE_FAILED", answer="A2 growth to large clusters by 1 s remains unanswered."
        ),
        Q6=dict(
            status="NOT_RUN_M0_GATE_FAILED",
            answer="A2 B2→A and B3→B branch-steering robustness remains unanswered; no full 6 s trajectories were launched.",
        ),
    )
    decision = dict(
        status="A2_RUNTIME_GATE_FAILED",
        finished_utc=datetime.now(timezone.utc).isoformat(),
        study="MODEL SENSITIVITY / ABLATION STUDY",
        A2_promoted=False,
        M0_gate=gate,
        questions=answers,
        dense_table=table,
        simulations_started=dict(
            A1_dense=2,
            A2_dense=2,
            A0=0,
            L0=0,
            primary_only=0,
            secondary_dt_half=0,
            formation_025s=0,
            growth_1s=0,
            full_6s=0,
        ),
        forbidden_inferences=[
            "Corrected or validated near-field force",
            "True SonoVue equilibrium spacing",
            "Branch control robust or sensitive based only on dense probes",
        ],
    )
    write(BASE / "final_decision/FINAL_DECISION.json", decision)
    md = [
        "# Stage 4M final decision",
        "",
        "**A2_RUNTIME_GATE_FAILED** — MODEL SENSITIVITY / ABLATION STUDY.",
        "",
        "The four requested 5 ms dense diagnostics and the existing-data comparison are evaluated below. No 0.25 s, 1 s or 6 s A2 campaign was launched. L0/A0 and primary-only trajectories were not rerun.",
        "",
        "| Case | Policy | Status | Steps | Wall time (s) | Step reduction | Wall speedup |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in table:
        md.append(
            f"| {row['case']} | {row['policy']} | {row['status']} | {row['accepted_steps']} | {row['wall_s']:.1f} | {row['step_reduction']:.3f}× | {row['wall_speedup']:.3f}× |"
        )
    md += [
        "",
        "Both A2 cases must individually achieve ≥5× by either measure and pass numerical/envelope checks. This is an engineering gate, not an accuracy threshold. Historical A0 used an earlier contact-solver implementation with the same physical resistance, constraints and tolerance; runtime ratios therefore include a backend-revision and machine-load caveat.",
        "",
    ]
    for key, row in answers.items():
        md += [f"## {key}: {row['status']}", "", row["answer"], ""]
        if key == "Q2":
            for case, r in row["cases"].items():
                md.append(
                    f"{case}: {r['step_reduction']:.3f}× fewer accepted steps, {r['wall_speedup']:.3f}× wall speedup."
                )
            md.append("")
    md += ["## Dense-envelope and force-removal observations", ""]
    for case in ("B2", "B3"):
        for policy in ("A1", "A2"):
            row = data[case][policy]
            f = force[case][policy]
            if "final_rg_m" in row:
                md.append(
                    f"{case}/{policy}: same-ID Rg {row['initial_rg_m'] * 1e6:.3f}→{row['final_rg_m'] * 1e6:.3f} µm; maximum Rg/initial {row['max_rg_ratio']:.4f}; endpoint rigid-fit residual {row['rigid_fit_residual_end_m'] * 1e6:.3f} µm. Effective/classical pair-norm ratio {f['initial']['effective_classical_ratio']:.4f}→{f['final']['effective_classical_ratio']:.4f}."
                )
    md += [
        "",
        "Native accepted-state graph turnover, censored edge lifetimes, unordered pair-distance quantiles, subfloor counts, contact QP statistics and projection diagnostics are retained under `comparisons/dense/` and each new run. Graphs and positions are never interpolated. Same-ID clusters were formed by legacy L0 evolution; preservation or rearrangement is not evidence of A2 formation.",
        "",
        "## Scope and reproducibility",
        "",
        "The frozen Stage 4L translation backend hash is `0b3e36eae784f064986caf371a201f1ded0be5e6970eef6f1f856c6dd89c9725`. The only production changes are secondary pair weighting and config forwarding. Each run contains source and input identities, all accepted states and force diagnostics. Two B3 startup attempts failed before integration because a shared compiled cache was read during compilation; each was restarted with an independent cache, using unchanged contact source/compiler flags. The generic smooth-removal banner archived by the first runner applies to A2; A1 is explicitly a hard subfloor switch, as clarified by its POLICY_SCOPE.json.",
        "",
        "Figures M1, M2, M3 and M5 are provided in PNG/PDF. M4, M6 and M7 are unavailable because their experiment gates did not open; no synthetic replacement data were created.",
        "",
        "[Doinikov & Bouakaz (2015)](https://doi.org/10.1103/PhysRevE.92.043001) and [Chen, Chen & Li (2020)](https://doi.org/10.1016/j.ultsonch.2019.104814) motivate unresolved close-separation and multi-bubble coupling. Neither paper validates a zero near-contact force or provides correction factors for this study.",
    ]
    text = "\n".join(md) + "\n"
    (BASE / "final_decision/FINAL_REPORT.md").write_text(text)
    doc = ROOT / "docs/STAGE4M_NEAR_CONTACT_SECONDARY_ABLATION.md"
    contents = doc.read_text().split("\n## Final measured outcome\n")[0]
    doc.write_text(
        contents + "\n## Final measured outcome\n\n" + text.replace("# Stage 4M final decision\n\n", "", 1)
    )
    print(json.dumps(dict(status=decision["status"], dense_table=table)), flush=True)


if __name__ == "__main__":
    main()

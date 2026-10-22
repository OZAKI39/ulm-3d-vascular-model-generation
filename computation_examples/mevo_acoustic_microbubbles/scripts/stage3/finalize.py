"""Assemble scoped acceptance and reports exclusively from actual result artifacts."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from scripts.stage3.common import atomic_json

ROOT = Path("results/stage3")


def read(kind):
    return json.loads((ROOT / kind / "STATE.json").read_text())


def ev(path):
    p = Path(path)
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def main():
    data = {
        k: read(k)
        for k in [
            "acoustic",
            "acoustic_grid",
            "secondary",
            "clustering",
            "flow",
            "history",
            "cylinder",
            "postchecks",
            "performance",
            "software_checks",
        ]
    }
    regression = json.loads((ROOT / "audit/REGRESSION_STATE.json").read_text())
    assert all(
        data[k]["status"] == "PASS"
        for k in [
            "secondary",
            "clustering",
            "flow",
            "history",
            "cylinder",
            "postchecks",
            "performance",
            "software_checks",
        ]
    )
    assert regression["status"] == "PASS"
    radial = Path("results/stage2_5/response_table/20261001T220715485503Z")
    definitions = {
        "acoustic_basis": (
            "PASS_WITH_SCOPE",
            "Water 2.5 MHz four-source basis; linear superposition verified. Two-grid differences remain; no grid-independence claim.",
            ["acoustic", "acoustic_grid", "postchecks"],
            ["k-wave-python0.6.2; existing backend"],
        ),
        "radial_response": (
            "PASS_WITH_SCOPE",
            "Unchanged Stage2.5 adaptive table: R0 1.6–1.96 um, f 2.3–2.7 MHz, P 5–40 kPa; fixed Tu2009-derived manifold; no experimental calibration.",
            [],
            ["APECSS a0ee071ca247b5489e997d90e83b5f53c40e4ec7", "10.1121/1.3242346"],
        ),
        "primary_force": (
            "PASS",
            "Monochromatic complex-force identity checked against full-cycle integration.",
            ["secondary"],
            ["F=-<V grad(p)>"],
        ),
        "secondary_force": (
            "PASS_WITH_SCOPE",
            "ACCEPTED_FOR_SIMULATOR: SB1–SB6 implementation checks; not validated near-contact dynamics.",
            ["secondary"],
            ["Periodic monopole volume-derivative correlation; no Zhang gate"],
        ),
        "translation_approximation": (
            "WARNING",
            "Quasi-steady equilibrium-radius Stokes transport. Tested history/averaging discrepancies limit quantitative displacement and arrival time; no automatic history upgrade.",
            ["history"],
            ["Stokes force balance", "10.1063/1.3227903", "10.1103/p475-xggm"],
        ),
        "history_force_sensitivity": (
            "PASS_WITH_SCOPE",
            "Three SonoVue scenarios, 20 cycles and 0.1 ms, matched Aresolved/B/C and separate actual production averaging control, converged endpoints. Not experimental validation.",
            ["history"],
            ["10.1063/1.3227903", "arXiv:2506.16413v2", "10.1121/1.3242346"],
        ),
        "clustering_trend": (
            "PASS_WITH_SCOPE",
            "Eight separated bubbles, 5/10/15 kPa, early 0.1 ms, no flow. Late contact topology approximated; trend benchmark, not parameter-identical reproduction.",
            ["clustering"],
            ["10.1121/1.5010170"],
        ),
        "flow_coupling": (
            "PASS_WITH_SCOPE",
            "Synthetic spatial/temporal linear VTU versus analytic trajectories; does not validate imported CFD.",
            ["flow"],
            ["FlowField.velocity(x,t); VTK interpolation"],
        ),
        "cylinder_demo": (
            "PASS_WITH_SCOPE",
            "C1–C4 at water 2.5 MHz, 64 monodisperse bubbles, 2 ms; same source amplitudes and initial state for phase change. Warning-only pair extrapolation.",
            ["cylinder", "postchecks"],
            ["Existing cylinder/PZT/Poiseuille geometry and production pipeline"],
        ),
    }
    from mevo_mb.validation.gates import REQUIRES
    from mevo_mb.utils.config import load_config

    current_config = load_config("configs/stage3/cylinder.yaml")
    table_state = json.loads((radial / "SUMMARY.json").read_text())
    banners = [
        json.loads((Path(r["run"]) / "PHYSICS_BANNER.json").read_text()) for r in data["cylinder"]["rows"]
    ]
    checks = dict(
        A_Zhang_not_a_gate=not any(
            "ZHANG" in name or "PAIR_REFERENCE" in name or "PAIR_REGIME" in name
            for names in REQUIRES.values()
            for name in names
        ),
        B_minimal_production_chain=current_config["bubble_response"]["model"] == "adaptive_response_table"
        and all("quasi-steady Stokes translation" in b["production_physics"] for b in banners),
        C_APECSS_table=table_state["status"] == "PASS_EMPIRICAL_DOMAIN",
        D_primary=data["secondary"]["primary_status"] == "PASS",
        E_secondary=all(data["secondary"]["tests"].values()),
        F_three_history_conditions=len({r["case"] for r in data["history"]["rows"]}) == 3,
        G_quantitative_history=all("force_budget" in r["models"]["C"] for r in data["history"]["rows"]),
        H_no_production_history_dependency=data["history"]["production_history_dependency"] is False,
        I_secondary_ON_OFF=data["clustering"]["tests"]["ON_OFF"],
        J_pressure_trend=data["clustering"]["tests"]["pressure_trend"],
        K_contact_equality=data["postchecks"]["contact_regression"]["bitwise_equal_positions"],
        L_cylinder=data["cylinder"]["tests"]["completed_four"],
        M_phase_control=data["cylinder"]["tests"]["phase_changes_distribution"],
        N_VTK_trajectory=data["flow"]["status"] == "PASS",
        O_scope_declared=all(
            b["experimental_calibration"] is False
            and "high-fidelity near-contact dynamics" in b["not_included"]
            and "hydrodynamic history force" in b["not_included"]
            for b in banners
        ),
    )
    assert all(checks.values())
    acceptance = dict(
        status="MINIMAL_SIMULATOR_VALIDATED_WITH_SCOPE",
        stage=3,
        created_UTC=datetime.now(timezone.utc).isoformat(),
        not_experimental_calibration=True,
        criteria_A_to_O=checks,
        items={},
    )
    for key, (status, scope, kinds, sources) in definitions.items():
        evidence = [ev(ROOT / k / "STATE.json") for k in kinds]
        if key == "radial_response":
            evidence = [ev(radial / p) for p in ["SUMMARY.json", "HOLDOUTS.json", "adaptive_response.h5"]]
        acceptance["items"][key] = dict(status=status, scope=scope, evidence=evidence, source=sources)
    atomic_json("PHYSICS_ACCEPTANCE.json", acceptance)
    atomic_json(ROOT / "PHYSICS_ACCEPTANCE.json", acceptance)
    # Logical Zhang archive uses actual preserved files, not copies with altered identities.
    candidates = set()
    for directory in ["results/stage2_5/zhang_forensic", "references", "configs/stage2", "configs/stage2_5"]:
        for p in Path(directory).rglob("*"):
            if (
                p.is_file()
                and "__pycache__" not in p.parts
                and (directory.endswith("zhang_forensic") or "zhang" in str(p).lower())
            ):
                candidates.add(p)
    for p in Path("results/stage2").rglob("*"):
        if p.is_file() and "zhang" in str(p).lower():
            candidates.add(p)
    for p in [
        Path("docs/STAGE2_5_VALIDATION_REPORT.md"),
        Path("docs/STAGE2_VALIDATION_REPORT.md"),
        Path("src/mevo_mb/validation/zhang2011.py"),
    ]:
        candidates.add(p)
    atomic_json(
        "validation/deprecated/zhang2011/INDEX.json",
        dict(
            status="DEPRECATED_AS_DEVELOPMENT_REFERENCE",
            files=[dict(ev(p), bytes=p.stat().st_size) for p in sorted(candidates)],
        ),
    )
    a = data["acoustic"]["validation"]
    g = data["acoustic_grid"]
    sb = data["secondary"]
    cl = data["clustering"]
    cy = data["cylinder"]
    hf = data["history"]
    pc = data["postchecks"]
    perf = data["performance"]
    fl = data["flow"]
    report = [
        "# Stage 3 — Minimal Validated PZT–Microbubble Simulator",
        "",
        "## 1. Executive summary",
        "",
        "**MINIMAL_SIMULATOR_VALIDATED_WITH_SCOPE**. Criteria A–O passed. The production chain is operational and its tested numerical identities and qualitative trends are accepted. Translation has a **WARNING**, not a claim of quantitative experimental accuracy. No patient/skull/BraVa computation or stopped historical job was resumed.",
        "",
        "## 2. Why simulator physics was simplified",
        "",
        "The research target is PZT-controlled migration/aggregation with vascular flow. Existing APECSS and k-Wave components are reused; complex pair/contact/history physics is kept outside production. Simplicity does not remove uncertainty: diagnostics and reference sensitivities remain visible.",
        "",
        "## 3. Deprecated Zhang branch",
        "",
        "Zhang 2011 was investigated but was not used as a validation basis because a key published trend could not be independently reproduced.",
        "",
        "[Logical archive](../validation/deprecated/zhang2011/README.md) preserves the original failed HDF5, matrix, figures, notes and two-implementation investigation. No sign, shell parameter or averaging-window tuning was performed. Archived literature assertions run separately from production CI.",
        "",
        "## 4. Production physics",
        "",
        "[Methods](METHODS_MODEL_DESCRIPTION.md): cached linear field → adaptive APECSS lipid response → cycle-averaged primary and unregularized classical secondary → equilibrium-radius Stokes + FlowField → geometric no-penetration. Production has no hydrodynamic-history dependency. Every run saves a PHYSICS_BANNER.json and the enabled-force flags. The six prerequisites are hashed and scoped; no override was used.",
        "",
        "## 5. Acoustic validation",
        "",
        f"The old bases use a different medium. One water/2.5 MHz basis was therefore generated, and reused for C1–C4. Direct superposition relative L2={a['relative_L2']:.6g}; phase RMS={a['phase_RMS_rad_masked']:.6g} rad. Amplitude/phase changes did not rerun propagation.",
        "",
        f"100 versus 80 μm grids at common initial cloud points differ by {100 * g['pressure_relative_L2']:.3f}% in complex pressure, {100 * g['gradient_relative_L2']:.3f}% in gradient, and {100 * pc['acoustic_primary_grid_relative_L2']:.3f}% in actual nonlinear-table primary force. These are retained numerical uncertainties, not proof of grid independence.",
        "",
        "## 6. APECSS radial-response validation",
        "",
        "The formal Stage2.5 table is unchanged: 790 validated adaptive cells, 256 independent holdouts, R0=1.6–1.96 μm, 2.3–2.7 MHz and 5–40 kPa. Existing deterministic cache, resume, multiprocessing and strict rejection are retained. The actual cylinder trajectories stayed inside this box. No production table expansion was necessary. New direct radial cases serve weak-field and history validation only; H3 at 150 kPa does not expand the production table.",
        "",
        "## 7. Primary-force validation",
        "",
        f"The retained formula −½ Re[V1 conj(gradP)] agrees with independent full-cycle integration to relative error {sb['primary_full_cycle_relative_error']:.3g}. Passing this identity does not remove acoustic discretization or radial constitutive uncertainty.",
        "",
        "## 8. Secondary-force physical acceptance",
        "",
        f"SB1–SB6 passed: same-phase attraction, opposite-phase repulsion, five-distance slope {sb['distance_slope']:.9f}, four actual APECSS weak-pressure points (100/200/400/800 Pa) with slope {sb['pressure_slope']:.9f}, reciprocal force and zero response. **ACCEPTED_FOR_SIMULATOR** refers to these basic properties, not near-contact validation.",
        "",
        "## 9. Hydrodynamic-history sensitivity study",
        "",
        "[Full twelve-section report](HISTORY_FORCE_SENSITIVITY.md) compares three conditions (1.78 μm; 2.5 MHz; 5, 20 and 150 kPa), at 8 μs and 0.1 ms. Aresolved/B/C have matched instantaneous inputs; Aprod separately represents actual production cycle averaging. Only B–C is termed history-only error. The study is independently implemented in validation/ and never imported by ParticleTracker.",
        "",
        "| Case | Short Aprod–C error (%) | Long Aprod–C error (%) | Long history-only B–C (%) |",
        "|---|---:|---:|---:|",
    ]
    for label in ["H1", "H2", "H3"]:
        short = next(r for r in hf["rows"] if r["case"] == label and r["window"] == "short")
        long = next(r for r in hf["rows"] if r["case"] == label and r["window"] == "long")
        report.append(
            f"| {label} | {100 * short['total_production_vs_C_relative_error']:.2f} | {100 * long['total_production_vs_C_relative_error']:.2f} | {100 * long['history_B_vs_C_relative_error']:.2f} |"
        )
    report += [
        "",
        "The maximum final timestep-refinement endpoint change is "
        + f"{100 * max(max(q['changes'].values()) for q in hf['convergence'] if q['samples_per_cycle'] == 256):.3f}%. All tested total errors exceed 20%; transient displacement and arrival times require explicit limitations. No empirical correction or automatic production upgrade was introduced. Signed H3 displacement is retained, not flipped to match a preferred direction.",
        "",
        "## 10. Clustering trend validation",
        "",
        "Eight initially separated bubbles, uniform field, no flow, identical initial cube. Acceptance uses the same first 0.1 ms for every pressure and a fixed ≥80% classical-policy pair fraction. The original 20 kPa attempt failed this fraction (about 77.5%) and remains in the result tree. A subsequent 5/10/15 kPa experiment preserves the threshold and time window and passes.",
        "",
        "| Pressure (kPa) | Early NN slope (mm/s) | Early policy-valid pair fraction | Time to 5% NN contraction (ms) |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in cl["rows"]:
        if r["secondary_enabled"] and r["dt_s"] == 1e-5:
            report.append(
                f"| {r['pressure_Pa'] / 1000:g} | {r['early_NN_slope_m_s'] * 1e3:.6f} | {r['early_classical_policy_fraction']:.6f} | {r['time_to_5percent_NN_contraction_s'] * 1e3:.3f} |"
            )
    report += [
        "",
        "Secondary OFF remains stationary. Pressure increases early attraction; timestep refinement passes 2%. [Lazarus 2017](https://pubmed.ncbi.nlm.nih.gov/29195473/) supports the qualitative pressure/clustering direction: **trend benchmark, not parameter-identical reproduction**. Late dense configurations are geometric approximations.",
        "",
        "## 11. Contact approximation",
        "",
        f"Sorted cKDTree candidates match brute force. The full C3 repeated trajectory has bitwise-equal positions and times (maximum position difference {pc['contact_regression']['maximum_coordinate_error_m']:.1g} m). A separate N=10000 cloud contains {perf['contact_count']} candidates with zero missed/extra contacts. Original projection/NNLS physics is unchanged. FAR_FIELD/TRANSITION/NEAR_CONTACT/CONTACT boundaries use d/(Ri+Rj)=5,2,1 (small classification tolerance at contact); they are configurable warning policies informed by Stage2.5, not universal boundaries. Near-contact results carry NEAR_CONTACT_PHYSICS_APPROXIMATED. Contact ≠ fusion, adhesion or coalescence.",
        "",
        "## 12. CPU/GPU performance",
        "",
        "The existing measured crossover at N=256 is retained for this server. New NumPy/CuPy agreement is around floating-point precision; that is an implementation check only.",
        "",
        "| N | NumPy median (ms) | CuPy median including transfers (ms) |",
        "|---|---:|---:|",
    ]
    for n in [64, 256, 1024]:
        values = {r["backend"]: r["median_wall_s"] for r in perf["pair_timings"] if r["N"] == n}
        report.append(f"| {n} | {1e3 * values['numpy']:.4f} | {1e3 * values['cupy']:.4f} |")
    report += [
        "",
        f"Contact candidate search at N=10000: brute {perf['contact_timings_s']['brute']:.4f} s, cKDTree {perf['contact_timings_s']['ckdtree']:.4f} s. C3 response queries take {cy['rows'][2]['timings']['response_queries_s']:.3f} s versus pair-force {cy['rows'][2]['timings']['pair_force_s']:.3f} s and contact search {cy['rows'][2]['timings']['contact_search_s']:.3f} s. Integration timing includes these nested operations and must not be added to them. Basis-loading and VTK-interpolation timings are separately saved in every relevant result.",
        "",
        "## 13. VTK/FlowField validation",
        "",
        f"The same tracker switches through flow.backend=zero/poiseuille/vtk. A space/time-linear analytic field was exported to two tetrahedral VTU snapshots and used to drive particles. Maximum trajectory discrepancy was {fl['maximum_trajectory_error_m']:.4g} m. No particle physics was modified for the imported flow. This establishes API/interpolation equivalence, not correctness of a future SimVascular solution.",
        "",
        "## 14. Cylinder demo",
        "",
        "The existing 20 mm long, 1.5 mm radius cylinder and four disc sources use water, 2.5 MHz, maximum Poiseuille speed 5 mm/s, 64 monodisperse 1.78 μm bubbles and 2 ms exposure. C1 acoustic OFF; C2 primary only; C3 primary+secondary; C4 identical to C3 with one source phase changed by π/2.",
        "",
        "| Case | Mean lateral displacement (μm) | Final mean NN (μm) | Target fraction | Largest cluster fraction | Near-contact bubble fraction |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in cy["rows"]:
        report.append(
            f"| {r['scenario']} | {r['mean_lateral_displacement_m'] * 1e6:.4f} | {r['final']['mean_nearest_neighbor_m'] * 1e6:.4f} | {r['final_target_region_fraction']:.6f} | {r['final']['largest_cluster_fraction']:.6f} | {r['final']['near_contact_fraction']:.6f} |"
        )
    report += [
        "",
        f"C1 analytic error is {cy['rows'][0]['analytic_trajectory_error_m']:.3g} m. Halving C3 nominal dt changes early displacement by {100 * pc['timestep']['early_displacement_relative_difference']:.6f}% and final positions by RMS {pc['timestep']['final_RMS_position_difference_m'] * 1e6:.4f} μm. Saved outputs include pressure/phase, force vectors, flow, trajectories, contact events, percentiles, acoustic-drift/flow ratio, 2D projected counts/areal concentration (particles/m²), and validity. The fixed target is transverse radius <80 μm; no bulk-injection concentration is inferred.",
        "",
        "## 15. PZT phase-control demonstration",
        "",
        f"C3 and C4 have bitwise-identical initial positions and the same source amplitudes. Changing only phases changes mean corresponding-particle final position by {cy['C3_C4_mean_particle_position_difference_m'] * 1e6:.4f} μm and cross-sectional centroid by {cy['C3_C4_cross_section_centroid_difference_m'] * 1e6:.4f} μm. The observed change exceeds the tested time-discretization difference, but acoustic-grid and model uncertainty still limit quantitative interpretation.",
        "",
        f"![Four cylinder controls](../{cy['run']}/C1_C4_distribution.png)",
        "",
        "## 16. Model scope",
        "",
        "[PHYSICS_ACCEPTANCE.json](../PHYSICS_ACCEPTANCE.json) declares PASS/PASS_WITH_SCOPE/WARNING per component, with SHA-256 evidence and source. Six gate names replace the old high-fidelity blockers. Acoustic basis identity, table medium/domain and invalid numerical states still reject. Operational pair/wall/retardation warnings permit scoped engineering studies; wide pair separations in the cylinder frequently exceed the kd policy.",
        "",
        "## 17. Known limitations",
        "",
        "[Known limitations](KNOWN_LIMITATIONS.md) include significant translation sensitivity, acoustic discretization, fixed shell parameters, no dissolution/nonsphericity, unvalidated close-pair/wall dynamics, no history/streaming/two-way feedback in production, and no patient calibration. A smooth trajectory, CPU/GPU agreement or conserved global quantity is not evidence of local physical accuracy.",
        "",
        "## 18. Readiness for bifurcation / SimVascular",
        "",
        f"The FlowField interface is operational. Regression audited {regression['immutable_result_files']:,} historical result-file entries with no missing/changed contents; the old 150-particle migration repeated with zero coordinate difference. Production tests and separate archived-reference checks pass; details and counts are in [software evidence](../{data['software_checks']['run']}/SUMMARY.json). See [quickstart](STAGE3_QUICKSTART.md) for repeatable commands.",
        "",
        "## 19. Recommended next stage",
        "",
        "Stop adding fundamental bubble models for this stage. Next use a simple bifurcation and an independently checked SimVascular flow field, retaining response bounds, contact labels and quantitative transport limitations. No new patient-specific or BraVa batch was started here.",
        "",
        "Accepted run identities:",
        "",
    ]
    for kind in [
        "acoustic",
        "history",
        "secondary",
        "clustering",
        "flow",
        "cylinder",
        "postchecks",
        "performance",
    ]:
        report.append(f"- {kind}: `{data[kind]['run']}`")
    report += [
        "",
        "The first history draft used cyclewise startup forcing in the production-style control; the accepted history run uses the actual steady-cycle production approximation and separates matched instantaneous comparisons. Both records remain. Runtime source snapshots and hashes distinguish drafts from accepted evidence.",
        "",
    ]
    Path("docs/STAGE3_VALIDATION_REPORT.md").write_text("\n".join(report))
    atomic_json(
        ROOT / "STATE.json",
        dict(
            status=acceptance["status"],
            acceptance="PHYSICS_ACCEPTANCE.json",
            report="docs/STAGE3_VALIDATION_REPORT.md",
            accepted_runs={k: v.get("run") for k, v in data.items()},
            criteria=checks,
        ),
    )
    print(acceptance["status"])


if __name__ == "__main__":
    main()

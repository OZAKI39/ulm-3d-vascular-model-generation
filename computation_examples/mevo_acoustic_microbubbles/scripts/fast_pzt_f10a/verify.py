"""Delivery checks on actual results; never reruns historical science."""

import subprocess
from datetime import datetime, timezone
import numpy as np
from .common import BASE, ROOT, js, write, sha


def main():
    stat = js(BASE / "static_v1/RESULTS.json")
    rad = js(BASE / "radial_v1/RESULTS.json")
    analysis = js(BASE / "analysis/RESULTS.json")
    rows = stat["rows"]
    for r in rows:
        assert r["status"] == "VALID_STATIC"
        assert len(r["per_bubble_pressure_Pa"]) == 32
        pp = np.array(r["per_bubble_pressure_Pa"])
        assert abs(np.median(pp) - r["pressure_Pa"]) < 1e-8
        assert 0.1 <= pp.min() <= pp.max() <= 40000
        assert max(r["event"]["amplitude_Pa"]) <= 20000 * (1 + 1e-8)
        if r["mode"] == "O":
            assert r["source_budget_Pa2"] <= r["S_budget_Pa2"] * (1 + 1e-8)
            assert pp.max() <= r["S_max_pressure_Pa"] * (1.01 + 1e-8)
    for name in ["middle", "interior"]:
        original = np.array(stat["states"]["near"]["positions_m"])
        translated = np.array(stat["states"][name]["positions_m"])
        np.testing.assert_allclose(translated - translated[0], original - original[0], rtol=0, atol=3e-19)
        assert stat["states"][name]["support_valid"]
    assert all(rad["pressure_admission"].values()) and rad["zero_exact"]
    assert all(v != "UNEXPECTED_ACCEPT" for v in rad["out_of_bounds"].values())
    assert all(r["RADIAL_NUMERICAL_PASS"] for r in rad["rows"])
    assert len(analysis["runs"]) == 13
    assert sum(r["quality"] == "PASS" for r in analysis["runs"].values()) == 11
    assert all(
        r["quality_raw"]["error_bits"] == 2048 for r in analysis["runs"].values() if r["quality"] == "FAIL"
    )
    assert analysis["flags"]["ACOUSTIC_EFFECT_RESOLVED"] and analysis["flags"]["CONTACT_DOMINATED"]
    assert not analysis["flags"]["TRAINING_READY"]
    for r in analysis["runs"].values():
        assert r["ledger_validation"]["status"] == "PASS"
    # All baseline tracked files are byte-identical; new F10A files only.
    baseline = "4bf369f4f65efb5909e9aedb69e460eb692584da"
    names = subprocess.check_output(
        ["git", "diff", "--name-only", baseline, "HEAD"], cwd=ROOT, text=True
    ).splitlines()

    def allowed(name):
        return (
            "/fast_pzt_f10a/" in name
            or "/fast_pzt/f10a/" in name
            or "/docs/FAST_PZT_F10A_" in name
            or "/docs/F10A_" in name
        )

    assert all(allowed(n) for n in names), names
    count = 0
    for name, rec in js(ROOT / "results/fast_pzt/f9/ARTIFACT_MANIFEST.json")["files"].items():
        assert sha(ROOT / "results/fast_pzt/f9" / name) == rec["sha256"]
        count += 1
    for name, rec in js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json").items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == rec["sha256"]
    for name, rec in js(BASE / "REMOTE_MANIFEST.json").items():
        assert sha(BASE / name) == rec["sha256"]
    for name, expected in js(BASE / "figures/FIGURE_SOURCES.json")["inputs"].items():
        assert sha(ROOT / name) == expected
    doc = js(BASE / "progress_document/UPDATE.json")
    assert doc["headings_unchanged"]
    assert sha(doc["path"]) == doc["after_sha256"]
    assert sha(doc["backup"]) == doc["before_sha256"]
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS"
    for name, expected in review["reviewed_sha256"].items():
        assert sha(BASE / "figures" / name) == expected
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            f9_manifest_files_verified=count,
            all_f5_inputs_verified=True,
            baseline=baseline,
            changed_paths=names,
            scope="New F10A paths only; old tracked files and all registered F9 scientific artifacts unchanged",
        ),
    )
    write(
        BASE / "analysis/ENGINEERING_CHECKS.json",
        dict(
            status="PASS",
            pressure_fairness_rows=len(rows),
            radial_cases=14,
            stored_attempts=13,
            quality_pass=11,
            expected_contact_fail=2,
            additional_infrastructure_failures=2,
            ledger_corruption_checks=analysis["ledger_corruption_checks"],
            geometry_rigid_translation=True,
            source_and_artifact_hashes=True,
            failures="pilot_v1 JSON serialization; pilot_v2 overly absolute arithmetic summation tolerance; corrected only in F10A. 3kPa residual-overlap failures deliberately not repaired by relaxing physics.",
            limitations="3kPa error2048 and rejected interval metadata retained; no accepted post-initial3kPa states. Rejected ledger vectors are zero/invalid and must not be used for a force-displacement estimate. The legacy fixed-target navigation fields in HDF summaries are not F10A admission criteria; analysis/RESULTS.json is authoritative for the incremental task.",
        ),
    )
    timing = dict(
        radial_process_s=rad["process_wall_s"],
        radial_solver_sum_s=rad["apecss_sum_solver_s"],
        static_process_s=stat["process_wall_s"],
        optimizer_s=stat["optimizer_s"],
        candidate_evaluations=stat["total_evaluations"],
        support_process_s=js(BASE / "support_v1/CHECKS.json")["process_wall_s"],
        trajectory_campaign_process_s=sum(
            js(BASE / q / "CAMPAIGN.json")["process_wall_s"]
            for q in ["pilot_v3", "extended_10_v1", "extended_5_v1"]
        ),
        trajectory_integrator_wall_s=sum(r["wall_s"] for r in analysis["runs"].values()),
        trajectory_output_s=sum(r["output_s"] for r in analysis["runs"].values()),
        plot_generation_s=js(BASE / "figures/PLOT_RUN.json")["wall_s"],
        analysis_s=analysis["analysis_wall_s"],
        meaning="Process and integrator wall times include host work; CUDA event times are not exclusive kernel time. Output is included in campaign process time. Literature/development/visual review not separately measured. Eleven complete trajectories are paired numerical conditions, not eleven independent biological samples.",
    )
    write(BASE / "WORKFLOW_TIMING.json", timing)
    write(
        BASE / "analysis/FINAL_GATE.json",
        dict(
            status="COMPLETE_WITH_SCIENTIFIC_LIMITATIONS",
            timestamp_UTC=datetime.now(timezone.utc).isoformat(),
            engineering="PASS",
            history_preservation="PASS",
            archive="PASS",
            figures="PASS",
            progress_document="PASS",
            lab_template="COMPLETE_NOT_PROVIDED",
            flags=analysis["flags"],
            scientific_scope="Resolved1kPa group positive-Z increment in one translated interior32-bubble state over10ms; per-particle packing not converged;3kPa stopped; no navigation/training authorization.",
            next="F10B investigate contact/per-particle timestep and mobility sensitivity at1kPa; obtain experimental parameters before stronger-field production navigation.",
            source_results_sha256=sha(BASE / "analysis/RESULTS.json"),
            figure_sources_sha256=sha(BASE / "figures/FIGURE_SOURCES.json"),
        ),
    )
    print("F10A delivery checks PASS; scientific limitations retained")


if __name__ == "__main__":
    main()

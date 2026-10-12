"""Delivery audit of actual execution, immutable history, figures and document edit."""

import subprocess
from pathlib import Path
import numpy as np
import pymupdf as fitz
from PIL import Image
from .common import ROOT, BASE, BASELINE, js, sha, write


def main():
    checks = {}
    original = js(BASE / "audit_v1/AUDIT.json")
    assert original["instrumentation_position_delta_m"] == 0
    assert original["original_saved_position_delta_m"] == 0
    assert original["result"]["reason"] == "LINE_SEARCH_LIMIT"
    oldaudit = js(BASE / "analysis/ORIGINAL_AUDIT.json")
    assert not oldaudit["linear"][0]["true_pass"] and not oldaudit["linear"][1]["true_pass"]
    assert all(
        r["eigen_min"] > 0 and r["preconditioner_eigen_min"] > 0 and r["nonpositive_curvature_count"] == 0
        for r in oldaudit["linear"]
    )
    finaltrials = [r for r in oldaudit["evaluations"] if r["trial"] and r["outer"] == 4]
    assert len(finaltrials) == 20 and all(
        r["stable_accept"] and not r["original_accept"] for r in finaltrials
    )
    checks["exact_original_failure_and_independent_diagnosis"] = "PASS"
    first = js(BASE / "corrected_v1/RESULTS.json")
    second = js(BASE / "optimized_v2/RESULTS.json")
    for result in [first, second]:
        assert result["numerical_admitted"] and len(result["rows"]) == 19
        assert sum(r["expected_quality_rejection"] for r in result["rows"]) == 2
        for r in result["rows"]:
            assert r["test_pass"]
            g = r["gpu"]
            m = r["metrics"]
            lin = np.array(g["linear_audit"])
            assert np.isfinite(lin).all() and np.all(lin[:, :2] <= lin[:, 2, None])
            assert g["conservative_equation_bound_m"] <= 1e-14
            assert m["max_position_vs_independent_m"] <= 1e-14 and m["max_force_vs_independent_N"] <= 1e-15
            assert m["max_equation_residual_m"] <= 1e-14 and m["net_force_N"] <= 1e-20
            assert m["active_sets_equal"] and m["edge_sets_equal"]
            if r["expected_quality_rejection"]:
                assert not g["converged"] and g["reason"] == "OVERLAP_BUDGET"
            else:
                assert (
                    g["converged"] and m["max_overlap_m"] <= 17.8e-9 and m["max_single_iteration_m"] <= 178e-9
                )
        assert len(result["invalid"]) == 3 and all(r["rejected"] for r in result["invalid"])
    for a, b in zip(first["rows"], second["rows"]):
        assert a["name"] == b["name"] and a["dt_s"] == b["dt_s"]
        for field in ["positions_m", "contact_force_N"]:
            np.testing.assert_array_equal(a["gpu"][field], b["gpu"][field])
    checks["references_original_budgets_and_optimization_bitwise_equivalence"] = "PASS"
    guards = js(BASE / "guards_v1/RESULTS.json")
    assert guards["status"] == "PASS" and guards["linear_failure_rejected"]
    assert len(guards["merit_rows"]) == 24 and all(r["decisions_equal"] for r in guards["merit_rows"])
    checks["independent_merit_and_incomplete_linear_failure_guard"] = "PASS"
    perf = js(BASE / "performance_v1/RESULTS.json")
    assert len(perf["samples"]) == 3
    for s in perf["samples"]:
        assert len(s["measurements"]) == 4
        for r in s["measurements"]:
            assert r["quality_pass"] and r["matched_full_wall_s"] >= r["solve_and_preparation_s"]
            assert r["solver_converged"] == (r["name"] != "original")
            if r["name"] == "cpu":
                assert r["components"]["gradient_bound_m"] <= 1e-14 and r["components"]["cpu_assembly_s"] > 0
    for name, reads in [("original", 35), ("corrected", 17), ("optimized", 12)]:
        assert perf["diagnostics"][name]["result"]["timing"]["scalar_reads"] == reads
    assert js(BASE / "analysis/PERFORMANCE_SUMMARY.json")["new_gpu_over_cpu_cost"] > 1
    checks["numerical_gate_precedes_three_matched_pairs_and_production_stop"] = "PASS"
    for p in BASE.glob("*_PROCESS.json"):
        d = js(p)
        assert d["exit_code"] == 0, p.name
        for r in d["resource_samples"]:
            assert len(r["processes"].splitlines()) <= 1, (p.name, r)
    assert not js(BASE / "ENVIRONMENT.json")["final_compute_processes"]
    for name, e in js(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert sha(BASE / name) == e["sha256"], name
    checks["all_remote_results_verified_and_GPU_not_contended"] = "PASS"
    prior = ROOT / "results/fast_pzt/f10i"
    for name, e in js(prior / "ARTIFACT_MANIFEST.json")["files"].items():
        assert sha(prior / name) == e["sha256"], name
    changes = subprocess.check_output(
        ["git", "diff", "--name-only", BASELINE, "--", "src", "scripts", "configs", "docs", "results"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    for name in changes:
        assert "f10j" in name.lower(), name
    checks["F1_F10I_tracked_history_unmodified"] = "PASS"
    screen = js(BASE / "analysis/TOPOLOGY_SCREEN.json")
    assert screen["new_trajectories"] == 0 and screen["persistence_s"] == 0.002
    b2 = [r for r in screen["rows"] if r.get("is_same_geometric_B2")]
    assert len(b2) == 3 and all(r["component_count_range"] == [1, 1] for r in b2)
    c = js(BASE / "analysis/TOPOLOGY_CANDIDATE.json")
    assert not c["B2_validation_executed"] and c["absolute_start_time_s"] == 0.4175
    with np.load(BASE / "analysis/HISTORICAL_F5_CANDIDATE.npz") as f:
        assert np.array_equal(f["particle_ids"], np.arange(32)) and f["active"].all()
    checks["historical_topology_scope_and_actual_candidate"] = "PASS"
    provenance = js(BASE / "figures/FIGURE_PROVENANCE.json")
    assert provenance["script_sha256"] == sha(ROOT / "scripts/fast_pzt_f10j/figures.py")
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert len(review["figures"]) == 4
    pdf_checks = []
    for fig in provenance["figures"]:
        for ext in ["png", "pdf"]:
            p = BASE / "figures" / (fig["figure"] + "." + ext)
            assert sha(p) == fig["files"][ext]
        for path, value in fig["source_sha256"].items():
            assert sha(BASE / path) == value, path
        png = BASE / "figures" / (fig["figure"] + ".png")
        with Image.open(png) as im:
            assert im.width >= 2000 and np.asarray(im).std() > 10
        assert review["figures"][fig["figure"]]["png_sha256"] == sha(png)
        pdf = BASE / "figures" / (fig["figure"] + ".pdf")
        assert pdf.read_bytes().startswith(b"%PDF")
        with fitz.open(pdf) as doc:
            pdf_text = "\n".join(page.get_text() for page in doc)
            fonts = doc.get_page_fonts(0)
            assert fonts and all(doc.extract_font(font[0])[3] for font in fonts)
        assert len(pdf_text) > 1000
        if fig["figure"].endswith("scope"):
            assert "PERSISTENT TOPOLOGY VALIDATION NOT RUN" in pdf_text
        pdf_checks.append(dict(file=pdf.name, text_characters=len(pdf_text), embedded_fonts=True))
    checks["four_actual_PNG_PDF_figures_provenance_and_visual_review"] = "PASS"
    progress = js(BASE / "progress_document/VERIFICATION.json")
    assert sha(Path(progress["path"])) == progress["after_sha256"]
    assert sha(Path(progress["backup"])) == progress["before_sha256"]
    assert progress["all_original_bytes_preserved"] and progress["heading_structure_unchanged"]
    checks["progress_backup_minimal_edit_and_headings"] = "PASS"
    gate = js(BASE / "analysis/FINAL_GATE.json")
    assert gate["GPU_B2_PRODUCTION_NOT_ADMITTED"] and not gate["training_ready"] and not gate["production"]
    assert (
        gate["contact_physics"] == "PHYSICS_UNCALIBRATED"
        and gate["persistent_topology"] == "PERSISTENT TOPOLOGY VALIDATION NOT RUN"
    )
    for name in ["DEVELOPMENT", "GPU_SOLVER_AUDIT", "PERFORMANCE"]:
        assert (ROOT / ("docs/FAST_PZT_F10J_" + name + ".md")).stat().st_size > 2000
    checks["layered_scientific_limits_and_reports"] = "PASS"
    write(
        BASE / "DELIVERY_VERIFICATION.json",
        dict(
            status="PASS",
            checks=checks,
            pdf_checks=pdf_checks,
            execution_commits_unchanged=gate["execution_commits"],
            script_sha256=sha(__file__),
        ),
    )
    files = {
        str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(BASE.rglob("*"))
        if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
    }
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(scope="All F10J deliverables except this self-referential manifest", files=files),
    )
    print("All delivery checks PASS;", len(files), "artifact hashes recorded.")


if __name__ == "__main__":
    main()

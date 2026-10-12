"""Read-only delivery verification; --seal also records checks and artifact hashes."""

import argparse
import json
import subprocess
import h5py
import numpy as np
import pymupdf
from PIL import Image
from .common import ROOT, BASE, BASELINE, js, sha, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--seal", action="store_true")
    args = ap.parse_args()
    checks = {}
    for name, entry in js(ROOT / "results/fast_pzt/f10j/ARTIFACT_MANIFEST.json")["files"].items():
        assert sha(ROOT / "results/fast_pzt/f10j" / name) == entry["sha256"], name
    for name, entry in js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json").items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == entry["sha256"], name
    checks["historical_artifacts_71_and_inputs_13"] = True
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", BASELINE, "--", "src"], cwd=ROOT, text=True
    ).splitlines()
    assert set(changed) == {
        prefix + "src/mevo_mb/fast_sim/f10i_gpu.py",
        prefix + "src/mevo_mb/fast_sim/f10k_scaling.py",
    }
    original = subprocess.check_output(
        ["git", "show", BASELINE + ":" + prefix + "src/mevo_mb/fast_sim/f10i_gpu.py"], cwd=ROOT, text=True
    )
    current = (ROOT / "src/mevo_mb/fast_sim/f10i_gpu.py").read_text()
    hook = "    # Original public research limit stays32. F10K's explicitly opted-in subclass\n    # raises only allocation scope; all equations, kernels and solver targets stay unchanged.\n    max_fixed_state_bubbles = 32\n\n"
    assert current.replace(hook, "").replace("n <= self.max_fixed_state_bubbles", "n <= 32") == original
    checks["default_model_unchanged_opt_in_allocation_hook_only"] = True
    for name, entry in js(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert sha(BASE / name) == entry["sha256"], name
    remote = js(BASE / "REMOTE_VERIFICATION.json")
    assert remote["status"] == "PASS" and len(remote["source_files_against_git_archive"]) == 5
    assert remote["final_compute_processes"] == ""
    checks["remote_77_files_five_source_commits_including_startup_failures"] = True
    recovery = js(BASE / "topology_recovery_v1/RESULTS.json")
    for name, value in recovery["preserved_files"].items():
        assert sha(BASE / "topology_dt10_v1" / name) == value
    assert js(BASE / "topology_dt10_v1_PROCESS.json")["exit_code"] == 1
    assert recovery["no_science_rerun"] and recovery["recovered_only_metadata"]
    checks["original_failed_exit_and_scientific_files_preserved"] = True
    analysis = BASE / "analysis"
    admission = js(analysis / "F5_ADMISSION.json")
    assert admission["admission"] == "B2_TOPOLOGY_CASE_NOT_ADMITTED"
    assert abs(admission["exact_wall_gap_min_m"]) < 1e-15
    alternatives = js(analysis / "ALTERNATIVE_SCREEN.json")
    selected = alternatives["selected"]
    assert sha(ROOT / selected["source"]) == selected["sha256"]
    assert selected["exact_wall_gap_min_m"] > 10 * selected["radius_m"]
    assert js(BASE / "topology_dt10_v1/RESTART_CHECK.json")["bitwise_positions"]
    t = js(analysis / "B2_TOPOLOGY_COMPARISON.json")
    c = t["comparison"]
    assert c["status"] == "MORPHOLOGY_NOT_CONVERGED" and not c["gates"]["edge_events"]
    assert c["gates"]["persistent_events"] and c["event_counts"] == [46, 47]
    assert sum(a["accepted_steps"] for a in t["audits"]) == 900
    for a in t["audits"]:
        assert a["contact_quality_pass"] and a["edge_ledger_exact"]
        assert a["max_equation_residual_m"] <= 1e-13
        assert a["max_overlap_m"] <= 17.8e-9 and a["max_single_correction_m"] <= 178e-9
        assert a["topology"]["counts"]["merge_persistent"] == 1
        assert a["topology"]["counts"]["split_persistent"] == 0
        assert sha(ROOT / a["path"] / "trajectory.h5") == a["trajectory_sha256"]
    assert js(analysis / "BOUNDARY_AUDIT.json")["status"] == "PASS"
    checks["admission_restart_contact_persistent_merge_and_retained_raw_edge_FAIL"] = True
    nav = js(analysis / "NAVIGATION.json")
    for row in nav["rows"]:
        assert sha(ROOT / row["source"]) == row["sha256"]
        with h5py.File(ROOT / row["source"]) as h:
            assert h["positions_m"].shape[1] == 1
            positions = h["positions_m"][:, 0]
            times = h["time_s"][:]
        assert positions.ndim == 3
        centroid = positions.mean(axis=1)
        cross = row["crossing"]
        k = cross["left_index"]
        assert centroid[k, 0] <= row["plane_x_m"] <= centroid[k + 1, 0]
        f = (row["plane_x_m"] - centroid[k, 0]) / (centroid[k + 1, 0] - centroid[k, 0])
        assert abs(times[k] + f * (times[k + 1] - times[k]) - cross["time_s"]) < 1e-14
        np.testing.assert_allclose(
            centroid[k] * (1 - f) + centroid[k + 1] * f, cross["centroid_m"], rtol=0, atol=1e-17
        )
        assert row["all_active"] and row["CFD_saved_valid"] and row["field_saved_valid"]
        assert all(row["proposed_target_support"]["inside"])
    assert not nav["navigation_established"] and not nav["original50nm_effect_gate_at_plane"]
    assert js(analysis / "NAVIGATION_PREREGISTRATION.json")["status"] == "PREREGISTERED_PROPOSAL_NOT_EXECUTED"
    checks["six_actual_trajectories_crossings_and_unexecuted_fixed_target_protocol"] = True
    s = js(analysis / "SCALING.json")
    assert s["solve_and_warmup_s"] <= 120 and s["stage_wall_s"] <= 300
    total = 0
    for row in s["rows"]:
        if row["N"] == 1024:
            assert row["status"] == "NOT_RUN" and row["reason"] == "STATIC_MEMORY_BUDGET"
            assert "CPU_over_GPU_speedup" not in row
            continue
        assert row["status"] == "PASS" and row["comparison"]["pass"]
        assert row["construction"]["default_limit_guard"]
        assert all(
            b["num_threads"] == 1 for b in row["construction"]["actual_BLAS"] if b["user_api"] == "blas"
        )
        for method in ["cpu", "gpu"]:
            ms = [m for m in row["measurements"] if m["method"] == method]
            assert [m["kind"] for m in ms].count("warmup") == 1
            assert [m["kind"] for m in ms].count("formal") == 2
            median = np.median([m["full_wall_s"] for m in ms if m["kind"] == "formal"])
            assert abs(median - row["formal_median_s"][method]) < 1e-15
            for m in ms:
                total += 1
                assert m["quality_pass"] and m["converged"] and m["full_wall_s"] < 10
                assert m["metrics"]["equation_residual_m"] <= 1e-14
                assert m["metrics"]["max_overlap_m"] <= 17.8e-9
                assert m["metrics"]["single_correction_m"] <= 178e-9
                if method == "gpu":
                    audit = np.array(m["components"]["linear_audit"])
                    assert np.all(audit[:, :2] <= audit[:, 2, None])
                else:
                    assert m["components"]["base_true_residual"] <= m["components"]["base_target"]
                    assert all(v["pass_"] for v in m["components"]["direction_linear_audit"])
        assert (
            abs(row["formal_median_s"]["cpu"] / row["formal_median_s"]["gpu"] - row["CPU_over_GPU_speedup"])
            < 1e-12
        )
    assert total == 18
    checks["18_same_quality_solves_threads_timing_budgets_and_1024_NOT_RUN"] = True
    provenance = js(BASE / "figures/FIGURE_PROVENANCE.json")
    for fig in provenance["figures"]:
        for name, value in fig["source_sha256"].items():
            assert sha(ROOT / name) == value, name
        for ext, value in fig["files"].items():
            path = BASE / "figures" / (fig["figure"] + "." + ext)
            assert sha(path) == value
            if ext == "png":
                with Image.open(path) as im:
                    im.verify()
            else:
                with pymupdf.open(path) as pdf:
                    assert len(pdf) == 1 and len(pdf[0].get_text()) > 500
                    assert all(font[2] == "Type0" for font in pdf[0].get_fonts())
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and review["actually_opened_PNG_and_rendered_PDF"]
    for name, value in review["files"].items():
        assert sha(BASE / "figures" / name) == value
    checks["four_PNG_PDF_sources_fonts_and_actual_visual_review"] = True
    gate = js(analysis / "FINAL_GATE.json")
    assert not gate["training_ready"] and gate["physics"] == "PHYSICS_UNCALIBRATED"
    assert gate["instantaneous_edge_gate"] == "FAIL" and gate["navigation"] == "NOT_ESTABLISHED"
    if not args.preflight:
        progress = js(BASE / "progress_document/CHECK.json")
        assert progress["status"] == "PASS" and progress["original_bytes_preserved"]
        for name, value in progress["file_sha256"].items():
            from pathlib import Path

            assert sha(Path(name)) == value
        checks["progress_backup_diff_headings_and_original_bytes"] = True
    if args.seal:
        assert not args.preflight
        write(BASE / "VERIFICATION.json", dict(status="PASS", checks=checks, identity=identity()))
        files = {
            str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(BASE.rglob("*"))
            if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
        }
        source_files = [
            *sorted((ROOT / "scripts/fast_pzt_f10k").glob("*.py")),
            *sorted((ROOT / "configs/fast_pzt/f10k").glob("*")),
            *sorted((ROOT / "docs").glob("FAST_PZT_F10K*.md")),
            ROOT / "src/mevo_mb/fast_sim/f10i_gpu.py",
            ROOT / "src/mevo_mb/fast_sim/f10k_scaling.py",
        ]
        write(
            BASE / "ARTIFACT_MANIFEST.json",
            dict(
                scope="F10K results; source and reports indexed separately; historical source files remain in baseline archives",
                files=files,
                source_sha256={str(p.relative_to(ROOT)): sha(p) for p in source_files if p.is_file()},
            ),
        )
    elif (BASE / "ARTIFACT_MANIFEST.json").exists():
        manifest = js(BASE / "ARTIFACT_MANIFEST.json")
        for name, v in manifest["files"].items():
            assert sha(BASE / name) == v["sha256"], name
        for name, v in manifest["source_sha256"].items():
            assert sha(ROOT / name) == v, name
    if args.preflight:
        write(BASE / "ENGINEERING_PREFLIGHT.json", dict(status="PASS", checks=checks, identity=identity()))
    print(json.dumps(dict(status="PASS", checks=checks), indent=2))


if __name__ == "__main__":
    main()

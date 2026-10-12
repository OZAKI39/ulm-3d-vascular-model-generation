"""Verify actual closed F9 delivery; never promote unresolved physical claims."""

from datetime import datetime, timezone
import difflib
import subprocess
import h5py
import numpy as np
from .common import ROOT, BASE, BASELINE, F5, js, sha, write
from .analyze import load_case


def verify(root, records):
    for name, row in records.items():
        path = root / name
        assert path.stat().st_size == row["bytes"] and sha(path) == row["sha256"], str(path)


def main():
    archive = js(BASE / "ARCHIVE_VERIFICATION.json")
    assert archive["status"] == "PASS"
    verify(BASE, archive["remote_files"])
    old_manifest = js(ROOT / "results/fast_pzt/f8/ARTIFACT_MANIFEST.json")
    verify(ROOT / "results/fast_pzt/f8", old_manifest["files"])
    inputs = js(F5 / "INPUT_MANIFEST.json")
    verify(F5, inputs)
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    inherited = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", BASELINE, "--", "src", "scripts", "configs"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    assert len(inherited) > 479
    for name in inherited:
        blob = subprocess.check_output(["git", "show", BASELINE + ":" + prefix + name], cwd=ROOT)
        assert (ROOT / name).read_bytes() == blob, name
    checks = js(BASE / "checks_v1/CHECKS.json")
    wall = js(BASE / "wall_admission_v1/CHECKS.json")
    assert checks["status"] == wall["status"] == "PASS"
    assert wall["audit_on_off_bitwise_equal"] and wall["F8_first_61_frames_bitwise_equal"]
    assert wall["actual_wall_event_count"] > 0 and wall["expected_wall_budget_failure_isolated"]
    assert js(BASE / "wall_admission_v1/EXPECTED_REJECTION.json")["status"] == "FAILED"
    with h5py.File(BASE / "wall_admission_v1/expected_rejection.h5") as h:
        assert len(h["rejected"]) == 1
    results = js(BASE / "analysis/RESULTS.json")
    assert set(results["cases"]) == {"Y0_10", "Y0_5", "Y1_10", "Y1_5", "Y2_10", "Y2_5"}
    count = events = 0
    for key in ["Y0_5", "Y1_10", "Y1_5", "Y2_10"]:
        folder = BASE / (key.lower() + "_v1")
        assert js(folder / "COMPLETE.json")["status"] == "COMPLETE"
        data = load_case(folder / "run")  # Recheck CPU closure against every accepted saved state.
        count += len(data["ledger"]["step_meta"])
        events += len(data["ledger"]["wall_contact_events"])
        reg = js(folder / "REGISTERED_PROTOCOL.json")
        historic = ROOT / reg["f8_source"]
        cfg, oldcfg = js(folder / "run/CONFIG.json"), js(historic / "CONFIG.json")
        assert cfg.pop("audit_enabled")
        if key == "Y0_5":
            assert cfg.pop("dt_s") == 5e-6 and oldcfg.pop("dt_s") == 1e-5
        assert cfg == oldcfg, key
        for name, expected in reg["source_hashes"].items():
            assert sha(historic / name) == expected
        assert (
            sha(ROOT / "results/fast_pzt/f8/y_static_v1/TARGETS_BEFORE_TRAJECTORIES.json")
            == reg["target_registration_sha256"]
        )
        if key != "Y0_5":
            reproduction = js(folder / "REPRODUCTION.json")
            assert all(
                reproduction[n]
                for n in ["positions_bitwise_equal", "gate_exact_equal", "actual_controls_equal"]
            )
            with h5py.File(folder / "run/trajectory.h5") as new, h5py.File(historic / "trajectory.h5") as old:
                np.testing.assert_array_equal(new["positions_m"][:], old["positions_m"][:])
                np.testing.assert_array_equal(new["time_s"][:], old["time_s"][:])
            assert js(folder / "run/GATE_CROSSING.json") == js(historic / "GATE_CROSSING.json")
            assert data["summary"]["controls"] == js(historic / "SUMMARY.json")["controls"]
    assert count == 180000 and events == 345231
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and len(review["figures"]) == 4
    for row in review["figures"]:
        stem = BASE / "figures" / row["name"]
        for ext in ["png", "pdf"]:
            assert sha(stem.with_suffix("." + ext)) == row[ext + "_sha256"]
        assert row["png_visually_opened"] and row["pdf_visually_opened_via_page_render"]
        assert not row["clipped_text_or_overlap"]
        provenance = js(stem.with_suffix(".json"))
        assert sha(stem.with_suffix(".json")) == row["provenance_sha256"]
        for name, expected in provenance["sources"].items():
            assert sha(ROOT / name) == expected, name
        for name, inputrow in provenance["input_hashes"].items():
            assert sha(F5 / name) == inputrow["sha256"]
        assert not review["pdf_checks"][row["name"]]["out_of_bounds"]
    update = js(BASE / "progress_document/UPDATE.json")
    assert update["status"] == "APPLIED" and update["headings_unchanged"]
    assert sha(update["path"]) == sha(BASE / "progress_document/after.md") == update["after_sha256"]
    assert sha(update["backup"]) == sha(BASE / "progress_document/before.md") == update["before_sha256"]
    before, after = [(BASE / "progress_document" / n).read_text() for n in ["before.md", "after.md"]]
    assert [s for s in before.splitlines() if s.startswith("#")] == [
        s for s in after.splitlines() if s.startswith("#")
    ]
    remaining = iter(after.splitlines())
    assert all(any(candidate == line for candidate in remaining) for line in before.splitlines())
    assert (BASE / "progress_document/UPDATE.diff").read_text() == "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile=update["backup"],
            tofile=update["path"],
        )
    )
    assert js(BASE / "analysis/AI_INTERFACE_CHECKS.json")["status"] == "PASS"
    ai = js(BASE / "analysis/AI_DATA_MANIFEST.json")
    assert all(not row["training_ready"] for row in ai["cases"].values())
    assert results["effects"]["Y0"]["TRANSVERSE_EFFECT_RESOLVED"]
    assert not results["effects"]["Y1"]["TRANSVERSE_EFFECT_RESOLVED"]
    assert js(BASE / "analysis/RPY_SCOPE.json")["classification"] == "NOT_ADMITTED"
    assert len(js(BASE / "JOB_TIMING.json")["jobs"]) == 6
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            baseline=BASELINE,
            inherited_source_files_unchanged=len(inherited),
            F8_files_rechecked=len(old_manifest["files"]),
            F5_inputs_rechecked=len(inputs),
            original_targets_budgets_and_controls_unchanged=True,
            historical_narrative_lines_preserved=True,
        ),
    )
    gate = js(BASE / "FINAL_GATE.json")
    assert gate["scientific"] == "WALL_SENSITIVE / NAVIGATION_NOT_RESOLVED"
    gate.update(
        figures="PASS",
        archive="PASS",
        progress_document="UPDATED",
        delivery_status="COMPLETE_WITH_RECORDED_LIMITATIONS",
        exact_ledger_accepted_steps=count,
        sparse_wall_events=events,
        new_numerical_controls=1,
        audit_reproductions=3,
        independent_new_initial_state_families=0,
        RPY_contact="NOT_ADMITTED",
        acceptance_scope="Software/ledger/noninterference and registered research deliverables complete. Geometry-wall sensitivity is not resolved physical accuracy. No navigation-success or AI-training approval.",
    )
    write(BASE / "FINAL_GATE.json", gate)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    sources = [
        *(ROOT / "scripts/fast_pzt_f9").glob("*.py"),
        *(ROOT / "src/mevo_mb/fast_sim").glob("f9_*.py"),
        ROOT / "configs/fast_pzt/f9/protocol.json",
        ROOT / "docs/FAST_PZT_F9_DEVELOPMENT.md",
        ROOT / "scripts/fast_pzt_f9/README.md",
    ]
    for p in sources:
        assert p.read_bytes() == subprocess.check_output(
            ["git", "show", commit + ":" + prefix + str(p.relative_to(ROOT))], cwd=ROOT
        ), p
    write(
        BASE / "DELIVERY_SOURCE.json",
        dict(
            source_commit=commit,
            source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources},
            scope="Committed delivery source; scientific execution commits remain in original IDENTITY files. Artifact commit follows this record.",
        ),
    )
    first = subprocess.check_output(
        ["git", "show", "-s", "--format=%cI", "687f8c33d"], cwd=ROOT, text=True
    ).strip()
    now = datetime.now(timezone.utc)
    write(
        BASE / "WORKFLOW_TIMING.json",
        dict(
            first_F9_commit_time=first,
            verification_time_UTC=now.isoformat(),
            recorded_workflow_span_s=(now - datetime.fromisoformat(first)).total_seconds(),
            scientific_process_sum_s=js(BASE / "JOB_TIMING.json")["total_process_elapsed_s"],
            scope="Elapsed workflow since first F9 commit; excludes initial reading/development before it, includes idle/analysis/transfer; not pure compute or total effort.",
        ),
    )
    manifest = dict(
        timestamp_UTC=now.isoformat(),
        scope="All F9 local scientific results, logs, analysis, figures and review/document evidence except this manifest; no dependency/cache builds.",
        files={
            str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(BASE.rglob("*"))
            if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
        },
    )
    write(BASE / "ARTIFACT_MANIFEST.json", manifest)
    verify(BASE, manifest["files"])
    print(
        dict(
            status=gate["delivery_status"],
            files=len(manifest["files"]),
            bytes=sum(row["bytes"] for row in manifest["files"].values()),
            inherited_source_files=len(inherited),
            ledger_steps=count,
            wall_events=events,
            figures=4,
        )
    )


if __name__ == "__main__":
    main()

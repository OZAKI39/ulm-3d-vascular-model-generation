"""Verify local delivery without running science or promoting scientific failures."""

from datetime import datetime, timezone
import difflib
import subprocess

from .common import BASE, BASELINE, F5, ROOT, js, sha, write


def verify_files(root, rows):
    for name, row in rows.items():
        p = root / name
        assert p.stat().st_size == row["bytes"] and sha(p) == row["sha256"], p


def main():
    archive = js(BASE / "ARCHIVE_VERIFICATION.json")
    assert archive["status"] == "PASS" and archive["remote_file_count"] == 100
    verify_files(BASE, archive["remote_files"])
    f7_manifest = js(ROOT / "results/fast_pzt/f7/ARTIFACT_MANIFEST.json")
    verify_files(ROOT / "results/fast_pzt/f7", f7_manifest["files"])
    inputs = js(F5 / "INPUT_MANIFEST.json")
    verify_files(F5, inputs)

    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and len(review["figures"]) == 4
    for fig in review["figures"]:
        stem = BASE / "figures" / fig["name"]
        for ext in ["png", "pdf"]:
            assert sha(stem.with_suffix("." + ext)) == fig[ext + "_sha256"]
        provenance = js(stem.with_suffix(".json"))
        assert sha(stem.with_suffix(".json")) == fig["provenance_sha256"]
        for name, expected in provenance["sources"].items():
            assert sha(ROOT / name) == expected, name
        for name, row in provenance["input_hashes"].items():
            assert sha(F5 / name) == row["sha256"], name

    update = js(BASE / "progress_document/UPDATE.json")
    assert update["status"] == "APPLIED" and update["headings_unchanged"]
    assert sha(update["path"]) == sha(BASE / "progress_document/after.md") == update["after_sha256"]
    assert sha(update["backup"]) == sha(BASE / "progress_document/before.md") == update["before_sha256"]
    before = (BASE / "progress_document/before.md").read_text()
    after = (BASE / "progress_document/after.md").read_text()

    def headings(s):
        return [line for line in s.splitlines() if line.startswith("#")]

    assert headings(before) == headings(after)
    assert (BASE / "progress_document/UPDATE.diff").read_text() == "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile=update["backup"],
            tofile=update["path"],
        )
    )
    # Apart from the documented F5–F7 -> F5–F8 input-range change, every historical line survives in order.
    old_lines = before.replace("F5–F7 实际复用的是 Stage 3", "F5–F8 实际复用的是 Stage 3", 1).splitlines()
    remaining = iter(after.splitlines())
    assert all(any(candidate == line for candidate in remaining) for line in old_lines)

    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", BASELINE, "--", "src", "scripts", "configs"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    assert len(paths) >= 479
    for name in paths:
        expected = subprocess.check_output(["git", "show", BASELINE + ":" + prefix + name], cwd=ROOT)
        assert (ROOT / name).read_bytes() == expected, name
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            baseline=BASELINE,
            unchanged_inherited_source_files=len(paths),
            F7_archived_files_rechecked=len(f7_manifest["files"]),
            original_inputs_rechecked=len(inputs),
            historical_narrative_lines_preserved=True,
        ),
    )

    results = js(BASE / "analysis/RESULTS.json")
    cases = results["cases"]
    assert set(cases) == {"A0", "A1", "Y0", "Y1", "Y2", "Y1_half", "Y2_half"}
    assert all(r["status"] == "COMPLETE" and r["numerical_quality"] == "PASS" for r in cases.values())
    assert all(r["all_IDs_active"] for r in cases.values())
    assert js(BASE / "analysis/AI_INTERFACE_CHECKS.json")["status"] == "PASS"
    ai = js(BASE / "analysis/AI_DATA_MANIFEST.json")
    for k, v in ai["cases"].items():
        assert v["navigation_label"] == cases[k]["final_evaluation"]["classification"]
        assert not v["training_ready"]
    assert ai["cases"]["Y0"]["navigation_label"] == "NAVIGATION_FAIL"
    gate = js(BASE / "FINAL_GATE.json")
    assert gate["software_checks"] == gate["A_only_factor_check"] == "PASS"
    assert gate["A1_bitwise_reproduces_F7"] and gate["A1_navigation_pass"] and not gate["A0_navigation_pass"]
    assert gate["Y_navigation"] == "COHESION_ONLY" and not gate["Y_resolved_additional_target_motion"]
    assert gate["contact"] == "NOT_ADMITTED" and gate["physics"] == "PHYSICS_UNCALIBRATED"
    gate.update(
        figures="PASS",
        archive="PASS",
        progress_document="UPDATED",
        delivery_status="COMPLETE_WITH_RECORDED_LIMITATIONS",
        new_completed_formal_trajectory_count=len(cases),
        new_completed_formal_trajectories=list(cases),
        new_completed_pilot_trajectory_count=5,
        all_produced_scientific_files_archived=True,
        acceptance_scope="Completed registered research, numerical checks, four reviewed PNG/PDF figures and WSL evidence. A1 no-flow passes its fixed-target task. Y target navigation not established; task-scale wall/step sensitivity unresolved. Contact cost not admitted. No physical validation or AI generalization claim.",
    )
    write(BASE / "FINAL_GATE.json", gate)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    write(
        BASE / "DELIVERY_SOURCE.json",
        dict(
            source_commit=commit,
            source_files={
                str(p.relative_to(ROOT)): sha(p)
                for p in sorted(
                    [
                        *(ROOT / "scripts/fast_pzt_f8").glob("*.py"),
                        *(ROOT / "src/mevo_mb/fast_sim").glob("f8_*.py"),
                        ROOT / "docs/FAST_PZT_F8_DEVELOPMENT.md",
                        ROOT / "scripts/fast_pzt_f8/README.md",
                        ROOT / "configs/fast_pzt/f8/protocol.json",
                    ]
                )
            },
            scope="Committed delivery source; archive commit is subsequent. Experiment commits remain in immutable IDENTITY records.",
        ),
    )
    first = subprocess.check_output(
        ["git", "show", "-s", "--format=%cI", "d00523a81"], cwd=ROOT, text=True
    ).strip()
    now = datetime.now(timezone.utc)
    write(
        BASE / "WORKFLOW_TIMING.json",
        dict(
            first_F8_commit_time=first,
            delivery_verification_time_UTC=now.isoformat(),
            recorded_workflow_span_s=(now - datetime.fromisoformat(first)).total_seconds(),
            scientific_process_sum_s=js(BASE / "JOB_TIMING.json")["total_process_elapsed_s"],
            scope="Workflow elapsed includes analysis, review, documentation and transfers after the first F8 commit; initial audit/development before that commit is not timed. Not active CPU/GPU cost or complete human effort.",
        ),
    )
    manifest = dict(
        timestamp_UTC=now.isoformat(),
        scope="All local F8 result artifacts except this manifest and dependency/cache builds; source and report are Git tracked.",
        files={
            str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(BASE.rglob("*"))
            if p.is_file()
            and p.name != "ARTIFACT_MANIFEST.json"
            and "dependencies" not in p.relative_to(BASE).parts
        },
    )
    write(BASE / "ARTIFACT_MANIFEST.json", manifest)
    verify_files(BASE, manifest["files"])
    print(
        dict(
            status=gate["delivery_status"],
            files=len(manifest["files"]),
            bytes=sum(v["bytes"] for v in manifest["files"].values()),
            inherited_files=len(paths),
            figures=4,
        )
    )


if __name__ == "__main__":
    main()

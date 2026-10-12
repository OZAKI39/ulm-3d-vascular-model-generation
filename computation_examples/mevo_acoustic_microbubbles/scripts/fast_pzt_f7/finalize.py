"""Close delivery gates after checking archived data and an existing visual review.

This does not run science or manufacture a visual-review pass. Regenerated plots
must first be inspected again and their review hashes updated deliberately.
"""

from datetime import datetime, timezone
import subprocess

from .common import BASE, BASELINE, F5, ROOT, js, sha, write


def main():
    archive = js(BASE / "ARCHIVE_VERIFICATION.json")
    assert archive["status"] == "PASS"
    for name, row in archive["remote_files"].items():
        p = BASE / name
        assert p.stat().st_size == row["bytes"] and sha(p) == row["sha256"], name

    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and len(review["figures"]) == 4
    for fig in review["figures"]:
        stem = BASE / "figures" / fig["name"]
        for extension in ["png", "pdf"]:
            assert sha(stem.with_suffix("." + extension)) == fig[extension + "_sha256"]
        provenance = js(stem.with_suffix(".json"))
        for name, expected in provenance["sources"].items():
            assert sha(ROOT / name) == expected, name
        for name, row in provenance["input_hashes"].items():
            assert sha(F5 / name) == row["sha256"], name

    update = js(BASE / "progress_document/UPDATE.json")
    assert update["status"] == "APPLIED" and update["headings_unchanged"]
    assert sha(update["path"]) == update["after_sha256"]
    assert sha(BASE / "progress_document/after.md") == update["after_sha256"]
    assert sha(update["backup"]) == update["before_sha256"]
    assert sha(BASE / "progress_document/before.md") == update["before_sha256"]

    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", BASELINE, "--", "src", "scripts", "configs"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    for name in paths:
        expected = subprocess.check_output(["git", "show", BASELINE + ":" + prefix + name], cwd=ROOT)
        assert (ROOT / name).read_bytes() == expected, name
    assert len(paths) == 479
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(status="PASS", baseline=BASELINE, unchanged_inherited_source_files=len(paths)),
    )

    results = js(BASE / "analysis/RESULTS.json")
    cases = results["cases"]
    completed = [name for name, row in cases.items() if row["status"] == "COMPLETE" and name != "old"]
    assert len(completed) == 7
    assert all(cases[name]["numerical_quality"] == "PASS" for name in completed)
    assert cases["slow"]["earlier_failed_prefix_bitwise_equal"]
    gate = js(BASE / "analysis/FINAL_GATE.json")
    for name in [
        "software_checks",
        "cached_predictor",
        "frozen_replay",
        "failure_retention_check",
        "timestep",
    ]:
        assert gate[name] == "PASS", name
    gate.update(
        figures="PASS",
        archive="PASS",
        progress_document="UPDATED",
        new_completed_full_trajectories=completed,
        new_completed_full_trajectory_count=len(completed),
        earlier_failed_attempts={
            "slow_v1": "Failed before HDF5 save; config, identity and error logs retained, original positions unavailable.",
            "slow_v2": "Saved 0–699.5 ms prefix and actual 700 ms failure state; not a completed trajectory.",
        },
        all_attempts_have_complete_trajectories=False,
        all_produced_scientific_files_archived=True,
        slow_full_run_status=cases["slow"]["status"],
        slow_numerical_quality=cases["slow"]["numerical_quality"],
        delivery_status="COMPLETE_WITH_RECORDED_LIMITATIONS",
        acceptance_scope=(
            "Software/numerical checks and seven completed F7 trajectories pass; "
            "earlier slow attempts remain failed, and slow_v1 has no saved trajectory. "
            "The completed slow_v3 passes numerical checks but fails navigation. "
            "No-flow v2 and its half-step check pass navigation. No physical validation."
        ),
    )
    write(BASE / "analysis/FINAL_GATE.json", gate)
    manifest = dict(
        timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        scope="All F7 result files excluding this manifest and dependency builds; source/report are Git tracked",
        files={
            str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(BASE.rglob("*"))
            if p.is_file()
            and p.name != "ARTIFACT_MANIFEST.json"
            and "dependencies" not in p.relative_to(BASE).parts
        },
    )
    write(BASE / "ARTIFACT_MANIFEST.json", manifest)
    for name, row in manifest["files"].items():
        p = BASE / name
        assert p.stat().st_size == row["bytes"] and sha(p) == row["sha256"]
    print(
        dict(
            status=gate["delivery_status"],
            files=len(manifest["files"]),
            bytes=sum(r["bytes"] for r in manifest["files"].values()),
            inherited_source_files=len(paths),
            figures=len(review["figures"]),
        )
    )


if __name__ == "__main__":
    main()

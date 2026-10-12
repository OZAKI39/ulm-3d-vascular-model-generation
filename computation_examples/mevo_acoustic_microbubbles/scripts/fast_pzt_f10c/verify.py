"""Traceable final audit; numerical PASS does not convert contested timing to PASS."""

from datetime import datetime, timezone
import hashlib
import json
import subprocess
import h5py
import numpy as np
from .common import ROOT, BASE, js, sha, write

BASELINE = "81fe1206d623f04ea7e3d3a78021251890185c9e"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"


def main():
    for stage in ["f10a", "f10b"]:
        folder = ROOT / "results/fast_pzt" / stage
        manifest = js(folder / "ARTIFACT_MANIFEST.json")["files"]
        for rel, meta in manifest.items():
            assert sha(folder / rel) == meta["sha256"], (stage, rel)
    for rel, meta in js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json").items():
        assert sha(ROOT / "results/fast_pzt/f5" / rel) == meta["sha256"]
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", BASELINE, "--", "src", "configs", "scripts", "docs", "results"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    unexpected = [name for name in changed if "f10c" not in name.lower()]
    assert not unexpected, unexpected
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            baseline_commit=BASELINE,
            changed_historical_files=[],
            historical_manifests_verified=["f10a", "f10b"],
            input_files_verified=13,
        ),
    )
    checked = set()
    for path in BASE.glob("*/IDENTITY.json"):
        identity = js(path)
        commit = identity["execution_commit"]
        for rel, digest in identity["source_sha256"].items():
            key = (commit, rel, digest)
            if key in checked:
                continue
            actual = subprocess.check_output(["git", "show", commit + ":" + PREFIX + rel], cwd=ROOT)
            assert hashlib.sha256(actual).hexdigest() == digest, key
            checked.add(key)
        registered = subprocess.check_output(
            ["git", "show", commit + ":" + PREFIX + "configs/fast_pzt/f10c/protocol.json"], cwd=ROOT
        )
        assert (
            hashlib.sha256(registered).hexdigest()
            == identity["protocol_sha256"]
            == sha(ROOT / "configs/fast_pzt/f10c/protocol.json")
        )
    write(
        BASE / "SOURCE_VERIFICATION.json",
        dict(
            status="PASS",
            unique_commit_file_hashes=len(checked),
            science_jobs=[p.parent.name for p in BASE.glob("*/IDENTITY.json")],
        ),
    )
    remote = js(BASE / "REMOTE_MANIFEST.json")
    for rel, meta in remote["files"].items():
        assert sha(BASE / rel) == meta["sha256"], rel
    assert js(BASE / "INTERFACE_CHECKS.json")["status"] == "PASS"
    assert js(BASE / "checks_v2/RESULTS.json")["status"] == "PASS"
    assert js(BASE / "failure_v1/COMPLETE.json")["accepted_prefix_bitwise"]
    result = js(BASE / "analysis/RESULTS.json")
    assert [r["B"] for r in result["batch"]] == [1, 4, 16, 64]
    assert all(
        r["status"] == "BATCH_NUMERICAL_PASS"
        and len(r["comparison"]) == r["B"]
        and all(c["bitwise"] for c in r["comparison"])
        for r in result["batch"]
    )
    assert result["decomposition"]["status"] == "PASS"
    assert all(not v["morphology_pass"] for k, v in result["modes"].items() if k != "Q0")
    assert js(BASE / "cold_v1/RESULTS.json")["status"] == "PASS"
    exclusive = result["exclusive_performance"]
    assert exclusive["completion"]["status"] == exclusive["resource_audit"]["status"] == "PASS"
    assert exclusive["completion"]["protocol_sha256"] == sha(
        ROOT / "configs/fast_pzt/f10c/exclusive_performance.json"
    )
    execution_commit = js(BASE / "exclusive_v1/IDENTITY.json")["execution_commit"]
    registered_exclusive = subprocess.check_output(
        ["git", "show", execution_commit + ":" + PREFIX + "configs/fast_pzt/f10c/exclusive_performance.json"],
        cwd=ROOT,
    )
    assert hashlib.sha256(registered_exclusive).hexdigest() == exclusive["completion"]["protocol_sha256"]
    samples = [
        json.loads(line) for line in (BASE / "exclusive_v1/RESOURCE_SAMPLES.jsonl").read_text().splitlines()
    ]
    assert len(samples) == exclusive["resource_audit"]["sample_count"]
    assert all(row["exclusive_sample"] and len(row["gpu_pids"]) <= 1 for row in samples)
    assert all(
        line.split()[0] == "mevo_fast_pzt_f10c_exclusive_v1"
        for row in samples
        for line in row["science_running"]
    )
    assert {row["phase"] for row in samples} >= {"batch", "cold"}
    assert not exclusive["resource_audit"]["before"]["gpu_pids"]
    assert not exclusive["resource_audit"]["after"]["gpu_pids"]
    assert [row["B"] for row in exclusive["batch"]] == [1, 4, 16, 64]
    assert [row["B"] for row in exclusive["cold"]] == [1, 4, 16, 64]
    assert all(
        row["status"] == "BATCH_NUMERICAL_PASS"
        and len(row["comparison"]) == row["B"]
        and all(c["bitwise"] for c in row["comparison"])
        for row in exclusive["batch"]
    )
    assert all(row["bitwise_all_samples"] for row in exclusive["cold"])
    speed_pass = all(row["resident_speedup"] >= 2 for row in exclusive["batch"] if row["B"] >= 16)
    assert speed_pass == exclusive["completion"]["engineering_2x_pass"]
    data = []
    for path in BASE.rglob("trajectory.h5"):
        with h5py.File(path) as h:
            meta = json.loads(h.attrs["metadata_json"])
            assert meta["training_ready"] is False and meta["split_group"] == "F10A_translated_original_Y32"
            assert meta["batch_size"] == h["positions_m"].shape[1] and h["positions_m"].shape[2:] == (32, 3)
            assert len(set(meta["sample_id"])) == meta["batch_size"]
            for group in [
                "source",
                "medium",
                "flow",
                "bubbles",
                "control",
                "trajectory",
                "morphology",
                "quality",
            ]:
                assert group in h
            assert np.isfinite(h["positions_m"][:]).all()
            data.append(
                dict(
                    path=str(path.relative_to(BASE)),
                    sha256=sha(path),
                    status=meta["status"],
                    valid_samples=meta["valid_sample_count"],
                    benchmark_only=True,
                    training_ready=False,
                    performance_scope=(
                        "EXCLUSIVE_GPU_SAMPLED"
                        if path.relative_to(BASE).parts[0] in ["exclusive_v1_batch", "exclusive_v1_cold"]
                        else "ORIGINAL_CAMPAIGN_SEE_RESOURCE_AUDIT"
                    ),
                )
            )
    write(
        BASE / "DATASET_ADMISSION.json",
        dict(
            datasets=data,
            independent_training_families_added=0,
            scope="Benchmark variants belong to the existing F10A family; no physical or morphology calibration claim.",
        ),
    )
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS"
    for rel, digest in review["files"].items():
        assert sha(BASE / "figures" / rel) == digest, rel
    figures = js(BASE / "figures/FIGURE_SOURCES.json")
    for rel, digest in figures["data_sha256"].items():
        assert sha(ROOT / rel) == digest, rel
    actual = subprocess.check_output(
        ["git", "show", figures["execution_commit"] + ":" + PREFIX + "scripts/fast_pzt_f10c/figures.py"],
        cwd=ROOT,
    )
    assert hashlib.sha256(actual).hexdigest() == figures["script_sha256"]
    progress = js(BASE / "progress_document/VERIFICATION.json")
    assert (
        sha(progress["backup"]) == progress["before_sha256"]
        and sha(progress["original_path"]) == progress["after_sha256"]
    )
    assert progress["headings_before"] == progress["headings_after"]
    supplement = js(BASE / "progress_document/exclusive_update/VERIFICATION.json")
    assert sha(supplement["backup"]) == supplement["before_sha256"]
    assert sha(progress["original_path"]) == supplement["after_sha256"]
    assert supplement["headings_before"] == supplement["headings_after"]
    gate = dict(
        status="F10C_DELIVERY_COMPLETE",
        timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        engineering_numerical_checks="PASS",
        batch_numerical="BATCH_NUMERICAL_PASS",
        batch_sizes=[1, 4, 16, 64],
        serial_batch_comparisons=85,
        comparison_tolerance="BITWISE_ZERO",
        decomposition="PASS",
        multirate={key: value["label"] for key, value in result["modes"].items()},
        multirate_production="NOT_ADMITTED",
        q4="NOT_RUN_MORPHOLOGY_GATE_FAILED",
        exclusive_2x_performance="THROUGHPUT_IMPROVED" if speed_pass else "THROUGHPUT_GATE_NOT_MET",
        exclusive_resource_audit="PASS_SAMPLED",
        exclusive_serial_batch_comparisons=85,
        exclusive_cold_sample_comparisons=85,
        exclusive_speed_ratios={str(r["B"]): r["resident_speedup"] for r in exclusive["batch"]},
        exclusive_resident_valid_samples_per_gpu_hour={
            str(r["B"]): r["valid_samples_per_gpu_hour"] for r in exclusive["batch"]
        },
        exclusive_cold_valid_samples_per_gpu_hour={
            str(r["B"]): r["valid_samples_per_gpu_hour"] for r in exclusive["cold"]
        },
        performance_scope="Full32 10ms benchmark, same-control family, single observation per size; no sustained-throughput or confidence-interval claim",
        shared_load_speed_ratios={str(r["B"]): r["resident_speedup"] for r in result["batch"]},
        shared_load_throughput_not_exclusive={
            str(r["B"]): r["valid_samples_per_gpu_hour"] for r in result["batch"]
        },
        physical_model_calibrated=False,
        particle_morphology_physically_validated=False,
        training_ready=False,
        source_identity="PASS",
        history_preservation="PASS",
        remote_hash_archive="PASS",
        figures="PASS",
        progress_document="PASS",
        remaining_required_performance_evidence=None,
        source_scope="Full32 individual DOFs; same-control-family batching; global fail-stop error semantics; all runtime original physics budgets retained.",
    )
    write(BASE / "analysis/FINAL_GATE.json", gate)
    files = {
        str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(BASE.rglob("*"))
        if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
    }
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            files=files,
            bytes=sum(v["bytes"] for v in files.values()),
            status="VERIFIED",
            archive_scope="All actual scientific outputs, failed attempts, provenance, necessary figures and reports; no new copies of historical large input fields.",
        ),
    )
    print(json.dumps(gate, indent=2))
    print("Archived files", len(files))


if __name__ == "__main__":
    main()

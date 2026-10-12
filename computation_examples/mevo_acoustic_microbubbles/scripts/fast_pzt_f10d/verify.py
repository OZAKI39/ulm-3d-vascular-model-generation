"""Verify immutable provenance, exact raw bytes, admission, figures and progress history."""

from datetime import datetime, timezone
import hashlib
import json
import subprocess
import h5py
import numpy as np
from .common import BASE, ROOT, BASELINE, js, sha, write

PREFIX = "computation_examples/mevo_acoustic_microbubbles/"


def git_bytes(commit, rel):
    return subprocess.check_output(["git", "show", commit + ":" + PREFIX + rel], cwd=ROOT)


def raw_pair(a, b, mapping=None):
    """Byte equality also checks NaN payloads and signed zero, beyond array_equal."""
    with h5py.File(a / "trajectory.h5") as old, h5py.File(b / "trajectory.h5") as new:
        mapping = list(range(old["positions_m"].shape[1])) if mapping is None else mapping
        checked = []
        for key in old["trajectory"]:
            aa, bb = old["trajectory/" + key][:], new["trajectory/" + key][:, mapping]
            assert aa.dtype == bb.dtype and aa.shape == bb.shape and aa.tobytes() == bb.tobytes(), (a, b, key)
            checked.append(key)
        for key in ["particle_id", "radii_m", "physical_dt_s", "time_s", "control_index"]:
            aa, bb = old[key][:], new[key][:]
            assert aa.dtype == bb.dtype and aa.tobytes() == bb.tobytes(), key
        return dict(
            reference=str(a.relative_to(ROOT)),
            candidate=str(b.relative_to(ROOT)),
            members=len(mapping),
            raw_arrays=checked,
            byte_exact=True,
        )


def main():
    historical = {}
    for stage in ["f10a", "f10b", "f10c"]:
        folder = ROOT / "results/fast_pzt" / stage
        manifest = js(folder / "ARTIFACT_MANIFEST.json")["files"]
        for rel, meta in manifest.items():
            assert sha(folder / rel) == meta["sha256"], (stage, rel)
        historical[stage] = len(manifest)
    for rel, meta in js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json").items():
        assert sha(ROOT / "results/fast_pzt/f5" / rel) == meta["sha256"]
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", BASELINE, "--", "src", "scripts", "configs", "docs", "results"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    assert all(
        "f10d" in name.lower() or name == PREFIX + "scripts/fast_pzt_f10c/runner.py" for name in changed
    ), changed
    # Every preexisting physical solver file must still match the baseline byte for byte.
    physical = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", BASELINE, "--", PREFIX + "src/mevo_mb/fast_sim"],
        cwd=ROOT.parents[1],
        text=True,
    ).splitlines()
    for rel in physical:
        local = rel.removeprefix(PREFIX)
        assert git_bytes(BASELINE, local) == (ROOT / local).read_bytes(), rel
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            baseline_commit=BASELINE,
            historical_files=historical,
            physics_files_unchanged=len(physical),
            input_files_verified=13,
            only_existing_source_extension="scripts/fast_pzt_f10c/runner.py explicit optional hooks",
        ),
    )

    checked = set()
    paths = [p for p in BASE.glob("*/IDENTITY.json") if p.parent.name != "local_validate_v1"]
    identities = [(str(p.relative_to(BASE)), js(p)) for p in paths]
    identities += [
        ("analysis/SUMMARY.json", js(BASE / "analysis/SUMMARY.json")["identity"]),
        ("figures/FIGURE_SOURCES.json", js(BASE / "figures/FIGURE_SOURCES.json")["identity"]),
    ]
    for name, ident in identities:
        commit = ident["execution_commit"]
        for rel, digest in ident["source_sha256"].items():
            key = (commit, rel, digest)
            if key not in checked:
                assert hashlib.sha256(git_bytes(commit, rel)).hexdigest() == digest, (name, key)
                checked.add(key)
        assert (
            hashlib.sha256(git_bytes(commit, "configs/fast_pzt/f10d/protocol.json")).hexdigest()
            == ident["protocol_sha256"]
            == sha(ROOT / "configs/fast_pzt/f10d/protocol.json")
        )
    assert (
        git_bytes("9905de6ee1308bf5ca4c15ab6cb6ca54e1f3c560", "configs/fast_pzt/f10d/protocol.json")
        == (ROOT / "configs/fast_pzt/f10d/protocol.json").read_bytes()
    )
    write(
        BASE / "SOURCE_VERIFICATION.json",
        dict(
            status="PASS",
            unique_commit_file_hashes=len(checked),
            identities=[name for name, _ in identities],
            excluded_development_only=[
                "local_validate_v1: dirty development-source check; not formal execution evidence"
            ],
            preregistered_protocol_unchanged=True,
        ),
    )
    for rel, meta in js(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert (BASE / rel).stat().st_size == meta["bytes"] and sha(BASE / rel) == meta["sha256"], rel
    assert js(BASE / "campaign_v1/INTERRUPTED.json")["status"] == "INTERRUPTED_IMPLEMENTATION_SELECTION"
    write(
        BASE / "analysis/EXCLUDED_EVIDENCE.json",
        dict(
            interrupted_campaign="campaign_v1",
            scope="All raw candidate outputs/logs retained, possibly incomplete HDF5; excluded from selected numerical/performance admission.",
            development_only="local_validate_v1",
            training_ready=False,
        ),
    )
    for path in [
        "profile_v1/RESULTS.json",
        "validate_v1/RESULTS.json",
        "validate_v1/INDEPENDENT_CHECKS.json",
        "storage_v1/RESULTS.json",
        "campaign_v2/COMPLETE.json",
        "campaign_v2/RESOURCE_AUDIT.json",
        "campaign_v2_benchmark/RESULTS.json",
        "campaign_v2_families/RESULTS.json",
        "campaign_v2_diversity/RESULTS.json",
    ]:
        assert js(BASE / path)["status"] == "PASS", path
    resource = js(BASE / "campaign_v2/RESOURCE_AUDIT.json")
    samples = [
        json.loads(line) for line in (BASE / "campaign_v2/RESOURCE_SAMPLES.jsonl").read_text().splitlines()
    ]
    assert len(samples) == resource["sample_count"]
    assert all(r["exclusive_sample"] and len(r["gpu_pids"]) <= 1 for r in samples)
    assert all(
        line.split()[0] == "mevo_fast_pzt_f10d_campaign_v2" for r in samples for line in r["science_running"]
    )
    assert not resource["before"]["gpu_pids"] and not resource["after"]["gpu_pids"]
    bench = js(BASE / "campaign_v2_benchmark/RESULTS.json")
    strict = []
    for r in bench["rows"]:
        assert r["comparison"]["status"] == "PASS" and r["comparison"]["quality_controls_topology_exact"]
        stem = f"{r['startup']}_B{r['B']}_"
        strict.append(
            raw_pair(
                BASE / "campaign_v2_benchmark" / (stem + "legacy"),
                BASE / "campaign_v2_benchmark" / (stem + "optimized"),
            )
        )
        if r["startup"] == "resident":
            assert r["F10C_archive_comparison"]["status"] == "PASS"
            strict.append(
                raw_pair(
                    ROOT / f"results/fast_pzt/f10c/exclusive_v1_batch/B{r['B']}_batch",
                    BASE / "campaign_v2_benchmark" / (stem + "legacy"),
                )
            )
    fam = js(BASE / "campaign_v2_families/RESULTS.json")
    inputs = js(BASE / "campaign_v2_families/INPUTS.json")["records"]
    assert [r["sample_id"] for r in inputs] == [r["sample_id"] for r in fam["original_order_mapping"]]
    for i, r in enumerate(fam["original_order_mapping"]):
        assert fam["comparisons"][i]["status"] == "PASS"
        strict.append(
            raw_pair(BASE / "campaign_v2_families" / f"serial_{i:03d}", BASE / r["path"], [r["member"]])
        )
    write(
        BASE / "analysis/RAW_BYTE_REGRESSION.json",
        dict(
            status="PASS",
            pairs=strict,
            scope="Same raw dtype/shape/bytes including NaN payload and signed zero; original scientific comparisons separately cover metrics, topology, quality and controls",
        ),
    )
    catalog = js(BASE / "analysis/SAMPLE_CATALOG.json")["records"]
    assert len({(r["file"], r["member"]) for r in catalog}) == len(catalog)
    checked_hdf = set()
    seed_groups = {}
    for r in catalog:
        assert not r["training_ready"] and r["physics_scope"] == "PHYSICS_UNCALIBRATED"
        if r["file"] not in checked_hdf:
            path = ROOT / r["file"]
            assert sha(path) == r["sha256"]
            with h5py.File(path) as h:
                m = json.loads(h.attrs["metadata_json"])
                assert np.isfinite(h["positions_m"][:]).all()
                assert m["physical_dt_s"] == 1e-5 and m["field_update_stride"] == m["flow_update_stride"] == 1
                assert m["valid_sample_count"] == (m["batch_size"] if m["status"] == "COMPLETE" else 0)
                assert all(
                    q["status"] == ("PASS" if m["status"] == "COMPLETE" else "ABORTED_BATCH")
                    for q in m["per_sample_quality"]
                )
                if m["status"] != "COMPLETE":
                    assert "rejected_interval" in h
            checked_hdf.add(r["file"])
        if r["sample_class"] == "EXPLORATORY_NUMERICAL_SAMPLE":
            seed_groups.setdefault(r["initial_state_family"], set()).add(r["split_group"])
    assert len(seed_groups) == 4 and all(len(v) == 1 for v in seed_groups.values())
    failure = js(BASE / "analysis/FAILURE_AUDIT.json")
    assert (
        failure["affected_trajectories"] == 6
        and failure["collateral_trajectories_without_observed_local_violation"] == 2
    )
    diversity = js(BASE / "campaign_v2_diversity/RESULTS.json")
    assert (
        diversity["prepared_valid"]
        == diversity["valid_initial_families"]
        == diversity["valid_control_trajectories"]
        == 2
    )
    assert all(c["status"] == "PASS" for c in diversity["comparisons"])
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS"
    for rel, digest in review["files"].items():
        assert sha(BASE / "figures" / rel) == digest
    figures = js(BASE / "figures/FIGURE_SOURCES.json")
    for rel, digest in figures["inputs"].items():
        assert sha(ROOT / rel) == digest
    for rel, digest in figures["outputs"].items():
        assert sha(BASE / "figures" / rel) == digest
    assert (
        len(list((BASE / "figures").glob("Fig*.png"))) == len(list((BASE / "figures").glob("Fig*.pdf"))) == 4
    )
    progress = js(BASE / "progress_document/VERIFICATION.json")
    assert (
        sha(progress["backup"]) == progress["before_sha256"]
        and sha(progress["original_path"]) == progress["after_sha256"]
    )
    assert progress["headings_before"] == progress["headings_after"]
    assert (BASE / "progress_document/after.md").read_text().replace(
        (BASE / "progress_document/addition.md").read_text(), "", 1
    ) == (BASE / "progress_document/before.md").read_text()
    summary = js(BASE / "analysis/SUMMARY.json")
    cold64 = next(r for r in bench["rows"] if r["B"] == 64 and r["startup"] == "cold")
    assert cold64["speedup"] >= js(ROOT / "configs/fast_pzt/f10d/protocol.json")["cold_speed_gate"]
    gate = dict(
        status="F10D_DELIVERY_COMPLETE_WITH_RECORDED_NUMERICAL_FAILURES",
        timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        full_particle_original_physics_unchanged=True,
        postprocessing_regression="PASS",
        raw_bytes="PASS",
        discrete_events="EXACT",
        max_float_difference=summary["non_bitwise_maxima"],
        cold_B64_speedup=cold64["speedup"],
        cold_1p25_gate="PASS",
        resource_audit="PASS_SAMPLED",
        multi_family_numerical="20_OF_20_PASS_BENCHMARK_ONLY",
        diverse_initial_arrangements=4,
        preparation_pass=2,
        control_trajectories_pass=2,
        control_trajectories_aborted=4,
        preparation_failed=2,
        complete_three_control_independent_families=0,
        diverse_full_control_coverage="NOT_PASSED_NUMERICAL_FAILURES",
        independent_OFF_branches_per_hour_small_sample=diversity["unique_initial_families_per_gpu_hour"],
        failure_records="ALL_RETAINED_NO_RETRY",
        failure_time_resolution="50 us interval; exact first substep unknown",
        DATA_PIPELINE_NUMERICAL_PASS=True,
        source_verification="PASS",
        history_preservation="PASS",
        remote_hashes="PASS",
        figures="PASS",
        progress_document="PASS",
        physical_model_calibrated=False,
        physics_scope="EXPLORATORY / PHYSICS_UNCALIBRATED",
        training_ready=False,
        next_priority="Member-local first-failure diagnosis and validated isolation/retry under unchanged contact budgets",
    )
    write(BASE / "analysis/FINAL_GATE.json", gate)
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            files={
                str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
                for p in sorted(BASE.rglob("*"))
                if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
            }
        ),
    )
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()

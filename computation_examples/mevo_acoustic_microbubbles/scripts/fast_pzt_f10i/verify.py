"""Delivery integration checks; no new physical solve or mutation of old evidence."""

import argparse
import subprocess
import numpy as np
import h5py
import pymupdf
from .common import BASE, ROOT, H, TRIALS, BASELINE, js, sha, write, read_steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seal", action="store_true", help="Write verification and final artifact hashes")
    args = ap.parse_args()
    audit = js(BASE / "BASELINE_AUDIT.json")
    assert sha(H / "ARTIFACT_MANIFEST.json") == audit["F10H_manifest_sha256"]
    history = js(H / "ARTIFACT_MANIFEST.json")["files"]
    for name, entry in history.items():
        assert sha(H / name) == entry["sha256"], name
    for name, entry in audit["inputs"].items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == entry["sha256"], name
    remote = js(BASE / "REMOTE_MANIFEST.json")["files"]
    for name, entry in remote.items():
        assert sha(BASE / name) == entry["sha256"], name
    assert js(BASE / "REMOTE_VERIFICATION.json")["status"] == "PASS"
    p = js(ROOT / "configs/fast_pzt/f10i/protocol.json")
    assert p["edge_time_tolerance_s"] == 1e-5 and p["persistence_s"] == 0.002
    assert p["kc_N_m"] == 0.04641319518378658
    assert p["gpu_comparison"]["overlap_m"] == 1.78e-8
    gate = js(BASE / "analysis/FINAL_GATE.json")
    events = js(BASE / "analysis/EVENT_MATCHES.json")
    assert gate["raw_edge_timing"] == "RAW_EDGE_TIMING_FAIL" and not gate["training_ready"]
    assert len(events["events"]) == 29 and events["timing_fail_count"] == 12
    assert len(events["cross_pair_inversions"]) == 1 and events["edge_id_sequence_pass"]
    for mode, path in TRIALS.items():
        t, x, records = read_steps(path)
        data = np.load(BASE / f"analysis/ACCEPTED_GEOMETRY_{mode}.npz")
        np.testing.assert_array_equal(t, data["time_s"])
        np.testing.assert_array_equal(x, data["positions_m"])
        np.testing.assert_array_equal(data["particle_id"], np.arange(32))
        assert all(s["accepted"] and not s["error_bits"] for s in records)
        gap = (
            np.linalg.norm(x[:, data["pairs"][:, 1]] - x[:, data["pairs"][:, 0]], axis=2) - 2 * p["radius_m"]
        )
        np.testing.assert_array_equal(gap, data["surface_gap_m"])
        np.testing.assert_array_equal(gap <= 0.1 * p["radius_m"], data["linked"])
        reconstructed = []
        for f, k in np.argwhere(data["linked"][1:] != data["linked"][:-1]):
            reconstructed.append(
                (round(float(t[f + 1]), 12), *map(int, data["pairs"][k]), int(data["linked"][f + 1, k]))
            )
        with h5py.File(path / "trajectory.h5") as h:
            stored = [
                (round(float(time), 12), int(i), int(j), int(on))
                for time, (_, i, j, on) in zip(
                    h["edge_events/time_s"][:], h["edge_events/batch_i_j_linked"][:]
                )
            ]
        assert sorted(reconstructed) == sorted(stored) and len(stored) == 29
    group = js(BASE / "analysis/GROUP_TOPOLOGY.json")
    assert all(v["members_equal"] for v in group["comparison"])
    for m in group["modes"].values():
        assert all(v == [list(range(32))] for v in m["members"])
        assert m["topology"]["events"] == []
        assert all(v == 0 for v in m["topology"]["counts"].values())
    assert gate["persistent_event_accuracy"] == "PERSISTENT_TOPOLOGY_NOT_EVALUATED"
    solver = js(BASE / "analysis/SOLVER_TOLERANCE.json")
    failures = []
    for row in solver["rows"]:
        assert row["reference"]["reference_pass"]
        for v in row["variants"]:
            result = v["result"]
            if not result["converged"]:
                failures.append((row["name"], v["tolerance_m"]))
            assert np.all(np.array(result["pair_force_N"]) >= 0)
            assert np.all(np.array(result["pair_force_N"])[np.array(result["geometric_gap_m"]) >= 0] == 0)
            assert v["vs_archived"]["active_contact_set_equal"] and v["vs_archived"]["edge_set_equal"]
    assert failures == [("DT10_0141", 1e-15)]
    prop = js(BASE / "analysis/TOLERANCE_PROPAGATION.json")
    path = ROOT / prop["trajectory"]["path"]
    t, x, records = read_steps(path)
    np.testing.assert_allclose(t, np.arange(201) * 1e-5, rtol=0, atol=1e-16)
    assert len(records) == 200 and all(q["accepted"] and not q["error_bits"] for q in records)
    hcfg = js(TRIALS["DT10"] / "CONFIG.json")
    cfg = js(path / "CONFIG.json")
    hcfg["duration_s"] = 0.002
    hcfg["protocol_sha256"] = sha(ROOT / "configs/fast_pzt/f10i/protocol.json")
    hcfg["options"]["numerical_mode"] = "B2_GEOM_K2_TOL14_DT10"
    assert cfg == hcfg, "Only duration, protocol identity and tighter numerical mode may differ"
    with h5py.File(path / "trajectory.h5") as h:
        assert not h.attrs["training_ready"]
        np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
        assert len(h["edge_events/time_s"]) == 23
    assert prop["same_dt_tolerance_comparison"]["matched_edge_max_time_difference_s"] == 0
    assert gate["finer_2p5us"] == "NOT_RUN"
    gpu = js(BASE / "gpu_v2/RESULTS.json")
    assert len(gpu["rows"]) == 19 and all(r["pass_"] for r in gpu["rows"][:-1])
    assert not gpu["rows"][-1]["gpu"]["converged"] and gpu["rows"][-1]["gpu"]["reason"] == "LINE_SEARCH_LIMIT"
    assert gpu["status"] == "GPU_B2_NOT_ADMITTED" and "performance_samples" not in gpu
    assert len([r for r in gpu["rows"] if r["contact_quality"] == "EXPECTED_BUDGET_REJECTION"]) == 2
    assert js(BASE / "gpu_v1/ENGINEERING_FAILURE.json")["completed_scientific_cases"] == 0
    fid = js(BASE / "figures/FIGURE_IDENTITY.json")
    assert len(fid["outputs"]) == 8
    for name, digest in fid["outputs"].items():
        assert sha(BASE / "figures" / name) == digest
    for name, digest in fid["inputs"].items():
        assert sha(BASE / name) == digest
    for pdf in (BASE / "figures").glob("Fig*.pdf"):
        with pymupdf.open(pdf) as doc:
            assert len(doc) == 1 and len(doc[0].get_text()) > 500
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and review["outputs"] == fid["outputs"]
    pd = BASE / "progress_document"
    meta = js(pd / "VERIFICATION.json")
    from pathlib import Path

    assert sha(Path(meta["path"])) == meta["after_sha256"] == sha(pd / "after.md")
    assert sha(Path(meta["backup"])) == meta["before_sha256"] == sha(pd / "before.md")
    assert meta["heading_structure_unchanged"] and meta["all_original_bytes_preserved"]
    changes = subprocess.check_output(
        ["git", "diff", "--name-status", BASELINE, "--", "."], cwd=ROOT, text=True
    )
    for line in changes.splitlines():
        status, name = line.split("\t")
        assert status == "A", line
        assert "f10i" in name.lower(), name
    for name in ["DEVELOPMENT", "EVENT_DIAGNOSTICS", "GPU_FEASIBILITY"]:
        assert (ROOT / f"docs/FAST_PZT_F10I_{name}.md").stat().st_size > 3000
    if args.seal:
        write(
            BASE / "DELIVERY_VERIFICATION.json",
            dict(
                engineering_delivery_status="PASS",
                scientific_global_pass=False,
                baseline_commit=BASELINE,
                verified_at_commit=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
                ).strip(),
                F10H_files_unchanged=len(history),
                original_inputs_unchanged=len(audit["inputs"]),
                remote_files_SHA_verified=len(remote),
                old_accepted_steps_independently_checked=1500,
                raw_events_independently_reconstructed=58,
                matched_events=29,
                new_accepted_steps=200,
                tighter_CPU_nonconvergence_preserved=failures,
                gpu_attempts=2,
                GPU_B2_NOT_ADMITTED=True,
                matched_cost_benchmark="NOT_RUN",
                figures_PNG_PDF=8,
                actual_visual_review="PASS",
                progress_backup_and_headings="PASS",
                original_tracked_files_modified=0,
                quality_and_timing_limits_unchanged=True,
                production_default_changed=False,
                training_ready=False,
            ),
        )
        files = {
            str(f.relative_to(BASE)): dict(bytes=f.stat().st_size, sha256=sha(f))
            for f in sorted(BASE.rglob("*"))
            if f.is_file() and f.name != "ARTIFACT_MANIFEST.json"
        }
        write(
            BASE / "ARTIFACT_MANIFEST.json",
            dict(
                scope="All F10I delivery files except this manifest itself",
                files=files,
                files_count=len(files),
                total_bytes=sum(v["bytes"] for v in files.values()),
            ),
        )
    elif (BASE / "ARTIFACT_MANIFEST.json").exists():
        for name, entry in js(BASE / "ARTIFACT_MANIFEST.json")["files"].items():
            assert sha(BASE / name) == entry["sha256"], name
    print(
        "PASS: archive integrity, all accepted clocks/IDs/edges, failed gates, independent reference evidence, figures and progress preservation."
    )


if __name__ == "__main__":
    main()

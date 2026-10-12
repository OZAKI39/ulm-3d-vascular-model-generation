"""Finite offline acceptance audit; never launches scientific integration."""

import json
import subprocess
import numpy as np
import h5py
from mevo_mb.fast_sim.f10e_isolation import isolate
from .common import BASE, ROOT, HISTORY, BASELINE, js, sha, write
from .prefix import compare


def unit_checks():
    records = [dict(sample_id=str(i), events=[dict(time_s=0)], original=i) for i in range(4)]
    calls = []

    def execute(group, parent):
        calls.append([r["original"] for r in group])
        return dict(
            status="FAILED" if any(r["original"] == 2 for r in group) else "COMPLETE",
            attempt_id=str(len(calls)),
        )

    result = isolate(records, "parent", execute)
    assert calls == [[0, 1], [2, 3], [2], [3]]
    assert [r["classification"] for r in result["leaves"]] == ["RECOVERED_BY_ISOLATION"] * 2 + [
        "FAILED_INDEPENDENT",
        "RECOVERED_BY_ISOLATION",
    ]
    unresolved = isolate(
        records[:2],
        "parent",
        lambda group, parent: dict(status="UNRESOLVED", attempt_id=group[0]["sample_id"]),
    )
    assert all(r["classification"] == "ABORTED_BATCH_UNRESOLVED" for r in unresolved["leaves"])
    for bad in [
        [],
        [records[0]],
        records[:2] + [records[0]],
        [records[0], dict(records[1], events=[dict(time_s=1)])],
    ]:
        try:
            isolate(bad, "parent", execute)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid isolation family accepted")
    return dict(
        status="PASS",
        routing=calls,
        unresolved_preserved=True,
        scope="Software stubs only: not new scientific trajectories",
    )


def main():
    out = BASE / "analysis"
    write(out / "ISOLATOR_CHECKS.json", unit_checks())
    manifest = js(HISTORY / "ARTIFACT_MANIFEST.json")
    assert all(sha(HISTORY / k) == v["sha256"] for k, v in manifest["files"].items())
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", BASELINE, "--", "src/mevo_mb/fast_sim"], cwd=ROOT, text=True
    ).splitlines()
    assert all(p.endswith(("f10e_contact.py", "f10e_isolation.py")) for p in changed), changed
    write(
        out / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            F10D_files_verified=len(manifest["files"]),
            original_fast_sim_files_modified=False,
            added_fast_sim_files=changed,
            runner_change="Optional registered dt/class/export keywords only; C0 regression validates inherited arithmetic",
        ),
    )
    regressions = []
    for p in sorted((BASE / "isolation_v1").glob("*/PREFIX_REGRESSION.json")):
        r = js(p)
        assert r["status"] == "PASS"
        regressions.append(str(p.relative_to(BASE)))
    for p in sorted((BASE / "contact_v1").glob("*/DIAGNOSTIC_PREFIX_REGRESSION.json")):
        r = js(p)
        assert r["status"] == "PASS"
        regressions.append(str(p.relative_to(BASE)))
    same = compare(
        BASE / "alternatives_v1/O_1000_1102_C1/trajectory.h5",
        BASE / "extension_v1/O_1000_1102_C1/trajectory.h5",
    )
    assert same["status"] == "PASS"
    write(out / "C1_DIAGNOSTIC_OFF_REGRESSION.json", same)
    assert len(regressions) == 6
    trees = js(BASE / "isolation_v1/RESULTS.json")["trees"]
    leaves = [r for t in trees for r in t["leaves"]]
    assert len(leaves) == 4 and all(r["classification"] == "FAILED_INDEPENDENT" for r in leaves)
    quality = []
    for p in sorted(BASE.glob("*_v1/*/trajectory.h5")):
        config = js(p.parent / "CONFIG.json")
        summary = js(p.parent / "SUMMARY.json")
        with h5py.File(p) as h:
            x = h["positions_m"][:, 0]
            dt = config["dt_s"]
            r = float(h["radii_m"][0])
            assert r == 1.78e-6 and np.array_equal(h["particle_id"][:], np.arange(32))
            np.testing.assert_array_equal(x[0], config["initial_positions_m"][0])
            assert np.all(np.isfinite(x)) and np.all(h["physical_dt_s"][:] == dt)
            d = np.linalg.norm(x[:, :, None] - x[:, None, :], axis=-1)
            d[:, np.arange(32), np.arange(32)] = np.inf
            assert np.min(d) >= 1.99 * r - 1e-18
            assert np.max(h["trajectory/diag"][:, :, :, 3]) <= 0.01 * r
            assert np.max(h["trajectory/diag"][:, :, :, 6]) <= 0.1 * r
            if summary["status"] == "COMPLETE":
                assert float(h["time_s"][-1]) == config["duration_s"]
                assert summary["valid_sample_count"] == 1
            else:
                assert "rejected_interval" in h and summary["valid_sample_count"] == 0
            quality.append(
                dict(
                    path=str(p.parent.relative_to(BASE)),
                    status=summary["status"],
                    sampled_max_overlap_nm=float(max(0, 2 * r - d.min()) * 1e9),
                    quality=summary["per_sample_quality"],
                )
            )
    assert len(quality) == 10
    write(out / "ACCEPTED_STATE_CHECKS.json", dict(status="PASS", rows=quality))
    resources = []
    for p in sorted(BASE.glob("*_PROCESS.json")):
        r = js(p)
        assert r["exit_code"] == 0
        assert all(len(s["processes"].splitlines()) <= 1 for s in r["resource_samples"])
        resources.append(
            dict(path=p.name, samples=len(r["resource_samples"]), status="AT_MOST_ONE_GPU_PID_OBSERVED")
        )
    assert js(BASE / "REMOTE_VERIFICATION.json")["status"] == "PASS"
    assert js(BASE / "ENVIRONMENT.json")["final_compute_processes"] == ""
    write(
        out / "RESOURCE_AUDIT.json",
        dict(
            status="NO_GPU_OVERLAP_OBSERVED",
            jobs=resources,
            isolation="Pre-launch/finish checks and integer Supervisor timestamps only; no continuous sampling",
            limits="0.5 s snapshots; GPU host PID names can be unavailable in the container. No CPU exclusivity or pure-kernel timing claim. Other CPU postprocessing jobs were left running.",
        ),
    )
    assert js(BASE / "figures/VISUAL_REVIEW.json")["status"] == "PASS"
    assert js(BASE / "progress_document/VERIFICATION.json")["status"] == "PASS"
    yield_data = js(out / "EFFECTIVE_YIELD.json")
    write(
        out / "FINAL_GATE.json",
        dict(
            status="F10E_DELIVERY_COMPLETE_WITH_RETAINED_NUMERICAL_AND_PHYSICAL_LIMITS",
            BATCH_ISOLATION_PASS=True,
            INDEPENDENT_NUMERICAL_PASS=dict(
                original_S_O_retries="0/4", historical_OFF="2/2", alternative_C1_O="1/1"
            ),
            CONTACT_DIAGNOSTIC_COMPLETE=True,
            ALTERNATIVE_NUMERICAL_SCHEME=dict(
                C1="10ms_PASS_ONE_CASE_EXPLORATORY",
                C2="5ms_WINDOW_ONLY",
                C3="5ms_WINDOW_ONLY",
                default_changed=False,
            ),
            SCIENTIFIC_TRIPLET_COMPLETE=False,
            complete_same_mode_triplets="0/4",
            PHYSICS_UNCALIBRATED=True,
            TRAINING_READY=False,
            training_ready=False,
            common_prefix_regressions=regressions + ["analysis/C1_DIAGNOSTIC_OFF_REGRESSION.json"],
            original_quality_budgets_preserved=True,
            original_fast_sim_unchanged=True,
            remote_hashes="PASS",
            history="PASS",
            figures="PASS",
            progress_document="PASS",
            yield_data=yield_data,
        ),
    )
    print(
        json.dumps(
            dict(
                status="PASS",
                scientific_runs=10,
                prefix_regressions=7,
                original_retries_pass=0,
                alternative_full_pass=1,
                complete_triplets=0,
            )
        )
    )


if __name__ == "__main__":
    main()

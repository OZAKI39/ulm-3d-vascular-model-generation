"""Offline scientific, source, visual and archive gates; no simulation reruns."""

import json
import subprocess
import h5py
import numpy as np
from pathlib import Path
from .common import BASE, HISTORY, ROOT, BASELINE, js, sha, write
from .analyze import gz


def main():
    protocol = js(ROOT / "configs/fast_pzt/f10f/protocol.json")
    gatechecks = js(BASE / "checks_v1/RESULTS.json")
    assert gatechecks["status"] == "PASS"
    assert js(BASE / "failure_checks_v1/RESULTS.json")["status"] == "PASS"
    rows = js(BASE / "short_v1/RESULTS.json")["rows"]
    assert len(rows) == 6
    checked = []
    for r in rows:
        path = BASE / "short_v1" / (r["case"] + "_" + r["mode"])
        config = js(path / "CONFIG.json")
        summary = js(path / "SUMMARY.json")
        ledger = gz(path / "SUBSTEP_LEDGER.json.gz")["steps"]
        trace = gz(path / "CONTACT_TRACE.json.gz")["traces"]
        dt = r["setting"]["dt_s"]
        iters = r["setting"]["contact_iterations"]
        assert len(ledger) == r["attempted_steps"] and len(trace) == min(32, len(ledger))
        assert all(t["accepted"] for t in ledger[:-1]) and not ledger[-1]["accepted"]
        assert all(t["error_bits"] == 0 for t in ledger[:-1]) and ledger[-1]["error_bits"] == 2048
        for i, t in enumerate(ledger):
            assert abs(t["start_s"] - i * dt) < 1e-14
            assert abs(t["attempted_end_s"] - (i + 1) * dt) < 1e-14
            assert len(t["round_residual_m"]) == iters
            assert max(t["round_max_correction_m"]) <= protocol["single_pass_budget_m"]
            assert t["quality"]["wall"]["max_pass_correction_m"] <= protocol["single_pass_budget_m"]
            if t["accepted"]:
                assert t["final_overlap_m"] <= protocol["residual_budget_m"]
                assert t["quality"]["max_accepted_penetration_m"] <= protocol["residual_budget_m"]
        assert ledger[-1]["final_overlap_m"] > protocol["residual_budget_m"]
        assert summary["field_update_stride"] == summary["flow_update_stride"] == 1
        assert summary["graph_chunk"] == 1 and not summary["training_ready"]
        assert summary["valid_sample_count"] == 0 and summary["status"] == "FAILED"
        assert summary["first_rejection_step"]["sample_id"] == summary["sample_id"][0]
        parent_hash = summary["initial_state_sources"][0]["parent_sha256"]
        selection = js(HISTORY / "SELECTION.json")
        chosen = next(v for v in selection.values() if v["case_id"] == r["case"])
        assert parent_hash == chosen["sha256"]
        parent = ROOT / chosen["relative_path"]
        np.testing.assert_array_equal(
            config["initial_positions_m"], js(parent / "CONFIG.json")["initial_positions_m"]
        )
        assert config["events"] == js(parent / "CONFIG.json")["events"]
        with h5py.File(path / "trajectory.h5") as h:
            x = h["positions_m"][:, 0]
            assert np.isfinite(x).all()
            d = np.linalg.norm(x[:, :, None] - x[:, None, :], axis=-1)
            d[:, np.arange(32), np.arange(32)] = np.inf
            assert d.min() >= 1.99 * protocol["radius_m"] - 1e-18
            assert np.all(h["trajectory/active"][:] == 1)
            np.testing.assert_array_equal(h["radii_m"][:], np.full(32, protocol["radius_m"]))
            np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
            np.testing.assert_array_equal(h["substeps/accepted_dt_s"][:], np.full(len(ledger) - 1, dt))
            assert h["last_legal_state"].attrs["time_s"] == ledger[-1]["start_s"]
            np.testing.assert_array_equal(
                h["last_legal_state/positions_m"][0], trace[-1]["initial_positions_m"]
            )
            np.testing.assert_array_equal(
                h["last_legal_state/positions_m"][:], h["rejected_interval/positions_m"][:]
            )
            np.testing.assert_array_equal(
                h["rejected_interval/proposal_m"][0], trace[-1]["proposal_positions_m"]
            )
            et = h["edge_events/time_s"][:]
            assert np.all(et <= ledger[-1]["start_s"] + 1e-14)
            assert np.all(np.abs(et / dt - np.round(et / dt)) < 1e-9)
        if r["mode"] == "B0":
            assert js(path / "B0_REGRESSION.json")["status"] == "PASS"
        checked.append(
            dict(
                case=r["case"],
                mode=r["mode"],
                accepted_substeps=len(ledger) - 1,
                rejected_proposal_sha256=sha(path / "trajectory.h5"),
                status="EVIDENCE_VALIDATED",
            )
        )
    old = js(BASE / "BASELINE_AUDIT.json")["historical_files_sha256"]
    for name, value in old.items():
        assert sha(ROOT / name) == value, name
    actual_e = Path(
        "/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10e_20261011/computation_examples/mevo_acoustic_microbubbles"
    )
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=actual_e, text=True).strip() == BASELINE
    assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=actual_e, text=True).strip()
    for name, value in old.items():
        assert sha(actual_e / name) == value, name
    # No preexisting engine source changed. Only the opt-in recorder is modified.
    changes = subprocess.check_output(
        ["git", "diff", "--name-only", BASELINE, "--", "src/mevo_mb/fast_sim"], cwd=ROOT, text=True
    ).splitlines()
    assert all(n.endswith("/f10f_smallsteps.py") for n in changes), changes
    verified = js(BASE / "REMOTE_VERIFICATION.json")
    assert verified["status"] == "PASS"
    for name, v in js(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert sha(BASE / name) == v["sha256"]
    inputs = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for name, v in inputs.items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == v["sha256"], name
    visual = js(BASE / "figures/VISUAL_REVIEW.json")
    assert visual["status"] == "PASS" and len(visual["files"]) == 8
    for name, v in visual["files"].items():
        assert sha(BASE / "figures" / name) == v["sha256"]
    progress = js(BASE / "progress_document/VERIFICATION.json")
    assert progress["original_bytes_preserved"] and progress["headings_preserved"]
    assert sha(progress["backup"]) == progress["before_sha256"]
    assert sha("/home/lzy/projects/microbubble_aggregation.md") == progress["after_sha256"]
    for p in BASE.glob("*_PROCESS.json"):
        d = js(p)
        assert d["exit_code"] == 0, p
        assert all(
            len([line for line in s["processes"].splitlines() if line.strip()]) <= 1
            for s in d["resource_samples"]
        )
    stop = js(BASE / "STOP_DECISION.json")
    assert not any(stop[k] for k in ["full_extension_admitted", "SS12_admitted", "batch_admitted"])
    result = dict(
        status="F10F_DELIVERY_COMPLETE_WITH_NEGATIVE_SMALLSTEPS_RESULT",
        engineering="PASS",
        original_defaults_and_history="PASS",
        strict_accepted_step_quality="PASS",
        original_failure_windows={"SS2": "0/2", "SS4": "0/2"},
        contact_admission="FAIL",
        group_centroid_prefix="PASS_PREPARATION_PREFIX_ONLY",
        paired_acoustic_delta="NOT_EVALUATED",
        morphology="MORPHOLOGY_NOT_CONVERGED",
        O_smallsteps_morphology="NOT_EVALUATED_NO_NONZERO_ACCEPTED_PREFIX",
        persistent_events="PARTIAL_OBSERVATION_RIGHT_CENSORED",
        performance="NOT_ADMITTED",
        SS12="NOT_ADMITTED",
        full_smallsteps="NOT_ADMITTED",
        batch="NOT_ADMITTED",
        numerical_complete_10ms_samples=0,
        complete_same_mode_triplets=0,
        PHYSICS_UNCALIBRATED=True,
        training_ready=False,
        EXPLORATORY=True,
        cases=checked,
        remote_hash_verification=verified,
        local_inputs_verified=len(inputs),
        historical_files_verified=len(old),
        figures="PASS",
        progress_document="PASS",
        effective_yield=js(BASE / "analysis/COST_AND_YIELD.json"),
    )
    write(BASE / "analysis/FINAL_GATE.json", result)
    files = {
        str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(BASE.rglob("*"))
        if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
    }
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            files=files,
            scope="Self-excluding local F10F delivery manifest; committed Git archive checked separately",
        ),
    )
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "status",
                    "engineering",
                    "contact_admission",
                    "morphology",
                    "performance",
                    "training_ready",
                ]
            },
            indent=2,
        )
    )
    print(
        "Verified",
        len(files),
        "local artifacts;",
        len(old),
        "protected history;",
        len(inputs),
        "local inputs",
    )


if __name__ == "__main__":
    main()

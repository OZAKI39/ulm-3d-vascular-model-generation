"""Read-only scientific/archive verification; failures must remain scientific failures."""

import argparse
import ast
import io
import re
import subprocess
import tarfile
import numpy as np
import h5py
from .common import ROOT, BASE, OLD, CONFIG, read, write, sha


def verify(final=False):
    checks = {}
    original = (
        ROOT.parents[1].parent / "mevo_fast_pzt_f10n_20261013/computation_examples/mevo_acoustic_microbubbles"
    )
    assert original.exists()
    assert (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=original, text=True).strip()
        == "e2dfb4f41cea33082583bbf8424a963ee194e85c"
    )
    assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=original, text=True).strip()
    manifest = read(original / "results/straight_tube_exp1/f10n/ARTIFACT_MANIFEST.json")
    for rel, record in manifest["files"].items():
        assert sha(original / rel) == record["sha256"], rel
    # Required protected source/data are bitwise identical in the new worktree.
    for rel in [
        "src/mevo_mb/fast_sim/f10h_contact.py",
        "src/mevo_mb/fast_sim/f10g_contact.py",
        "src/mevo_mb/fast_sim/pdms_cavity_engine.py",
        "src/mevo_mb/coupling/pdms_cavity.py",
        "src/mevo_mb/acoustics/water_disc_proxy.py",
        "src/mevo_mb/fast_sim/cubic_field.py",
    ]:
        assert sha(ROOT / rel) == sha(original / rel), rel
    assert ast.literal_eval((BASE / "cleanup/F10N_VERIFY.json").read_text())["status"] == "PASS"
    subprocess.run(["git", "merge-base", "--is-ancestor", "19c11e885", "3f5a117c1"], cwd=ROOT, check=True)
    checks["cleanup_precedes_bridge_and_original_F10N_unchanged"] = True
    reg = read(BASE / "PREREGISTRATION.json")
    p = read(CONFIG)
    assert reg["protocol"] == p and reg["target_locked_before_new_trajectory"]
    assert p["XY_targets_m"] == {"X": [2e-7, 0, 0], "Y": [0, 2e-7, 0]}
    assert p["kc_N_m"] == 0.04641319518378658 and p["CPU_contact_atol_m"] == 1e-13
    for rel, digest in reg["input_sha256"].items():
        assert sha(OLD / rel) == digest, rel
    checks["original_inputs_full_precision_controls_and_preregistered_targets"] = True
    for rel, record in read(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert sha(BASE / rel) == record["sha256"], rel
    checks["all_62_remote_artifacts_exact_hashes"] = True
    total = 0
    failed = []
    count = 0
    states = np.load(OLD / "initial_states.npz")["aggregate_m"]
    for attempt in ["pilot_v1", "validation_v1"]:
        camp = read(BASE / attempt / "CAMPAIGN.json")
        ident = camp["identity"]
        total += camp["wall_s"]
        deployment = read(BASE / "execution" / f"{attempt}.json")
        assert deployment["sha"] == ident["execution_commit"]
        package = BASE / "execution" / (deployment["sha"] + ".tar.gz")
        assert sha(package) == deployment["archive_sha256"]
        with tarfile.open(fileobj=io.BytesIO(package.read_bytes())) as arc:
            for rel, digest in ident["source_sha256"].items():
                import hashlib

                assert hashlib.sha256(arc.extractfile(rel).read()).hexdigest() == digest, rel
        for r in camp["runs"]:
            count += 1
            folder = BASE / attempt / "runs" / r["name"]
            inp = read(folder / "INPUT.json")
            assert inp["preregistration_sha256"] == sha(BASE / "PREREGISTRATION.json")
            e = inp["control"]
            original_event = read(OLD / "campaign_v1/runs" / inp["source_F10N"] / "INPUT.json")["events"][0]
            for key in ["amplitude_Pa", "phase_rad", "enabled", "flow_enabled"]:
                assert e[key] == original_event[key]
            assert not e["flow_enabled"]
            assert sha(folder / "trajectory.h5") == r["trajectory_sha256"]
            with h5py.File(folder / "trajectory.h5") as h:
                x = h["positions_m"][:]
                t = h["time_s"][:]
                m = h["step_metrics"][:]
                assert np.isfinite(x).all() and np.all(np.diff(t) > 0)
                np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
                np.testing.assert_array_equal(x[0], states)
                assert len(x) == r["accepted_steps"] + 1 and len(m) == r["attempted_steps"]
                assert np.all(abs(x[:, :, 0]) + p["radius_m"] < 0.01625)
                assert np.all(np.linalg.norm(x[:, :, 1:], axis=-1) + p["radius_m"] < 0.0008)
                assert h.attrs["flow"] == "STRICT_ZERO"
                pp = abs(h["force_samples/pressure_Pa"][:])
                assert np.all((pp == 0) | ((pp >= 0.1) & (pp <= 40000)))
                accepted = m[:, 0].astype(bool)
                assert np.all(m[accepted, 1] <= 1.01e-13)
                assert np.all(m[accepted, 2] <= 0.01 * p["radius_m"])
                assert np.all(m[accepted, 3] <= 0.1 * p["radius_m"]) and np.count_nonzero(m[:, 5]) == 0
                assert np.isclose(t[-1], r["completed_s"], rtol=0, atol=1e-15)
                if r["status"] == "PASS":
                    assert accepted.all() and np.isclose(t[-1], inp["duration_s"], rtol=0, atol=1e-15)
                else:
                    failed.append(r["name"])
                    assert r["reasons"] == ["RESIDUAL_OVERLAP_BUDGET"]
                    assert not accepted[-1] and accepted[:-1].all()
                    proposal = h["failure_proposal_m"][:]
                    assert not np.array_equal(proposal, x[-1]) and m[-1, 2] > 0.01 * p["radius_m"]
                if inp["source_F10N"] == "OFF":
                    np.testing.assert_array_equal(x, np.broadcast_to(x[0], x.shape))
    assert count == 12 and failed == ["STRONG_DT1", "STRONG_DT05", "STRONG_FINE"] and total < 300
    checks["immutable_executions_12_runs_original_budgets_and_preserved_failures"] = True
    g = read(BASE / "analysis/FINAL_GATE.json")
    com = read(BASE / "analysis/COMPARISONS.json")
    assert not g["strong_B2_2ms_quality_pass"] and not g["strong_B2_2ms_fixed_id_pass"]
    assert g["strong_B2_mathematical_pass"] and g["weak_recheck_status"] == "PASS"
    assert not g["X_target_reached"] and not g["Y_target_reached"] and not g["weak_target_reached"]
    assert not g["training_ready"] and not g["physical_calibration_pass"]
    assert all(g["windows"][k] == "NOT_RUN" for k in ["5ms", "10ms", "30ms"])
    assert all(not com["STRONG_" + k]["fixed_id_pass"] for k in ["DT05", "FINE"])
    assert all(com["WEAK_" + k]["fixed_id_pass"] for k in ["DT05", "FINE"])
    assert com["WEAK_M0_VS_B2"]["fixed_id_pass"] is None
    checks["scientific_failure_not_promoted_to_navigation_or_AI_pass"] = True
    if final:
        prov = read(BASE / "figures/FIGURE_PROVENANCE.json")
        review = read(BASE / "figures/VISUAL_REVIEW.json")
        assert len(prov["figures"]) == 4 and review["status"] == "PASS"
        assert sha(ROOT / "scripts/straight_tube_exp1/f10o/figures.py") == prov["source_script_sha256"]
        for rel, digest in prov["input_sha256"].items():
            assert sha(ROOT / rel) == digest, rel
        for rec in prov["figures"]:
            for ext, digest in rec["file_sha256"].items():
                assert sha(BASE / "figures" / (rec["figure"] + "." + ext)) == digest
        for r in review["figures"]:
            assert r["opened_png"] and r["opened_pdf_page_preview"]
            for ext, digest in r["file_sha256"].items():
                assert sha(BASE / "figures" / (r["figure"] + "." + ext)) == digest
        checks["four_opened_data_driven_PNG_and_PDF_pairs"] = True
        d = BASE / "progress_document"
        a = (d / "before.md").read_bytes()
        b = (d / "after.md").read_bytes()
        r = read(d / "CHECK.json")
        assert sha(r["path"]) == r["after_sha256"] and sha(r["backup"]) == r["before_sha256"]
        i = b.index("\n**F10O ".encode())
        j = b.index("\n---\n\n### 尾声".encode(), i)
        assert b[:i] + b[j:] == a
        assert re.findall(rb"^#{1,6} .*$", a, re.M) == re.findall(rb"^#{1,6} .*$", b, re.M)
        checks["progress_original_bytes_formulas_and_heading_structure_preserved"] = True
    return dict(
        status="PASS",
        checks=checks,
        scope="Engineering and evidence integrity only; strong scientific admission remains FAIL",
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--final", action="store_true")
    a = p.parse_args()
    result = verify(a.final)
    write(
        BASE / "analysis" / ("FINAL_VERIFICATION.json" if a.final else "ENGINEERING_VERIFICATION.json"),
        result,
    )
    print(result)


if __name__ == "__main__":
    main()

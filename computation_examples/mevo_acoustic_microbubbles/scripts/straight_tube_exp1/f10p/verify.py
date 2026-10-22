"""Read-only checks of the scientific archive; admission is not upgraded by this check."""

import argparse
import copy
import hashlib
import re
import subprocess
import tarfile
import numpy as np
import h5py
from scipy.spatial.distance import pdist
from .common import ROOT, BASE, OLD, CONFIG, read, write, sha
from .analyze import layered


def verify(final=False):
    checks = {}
    p = read(CONFIG)
    original = ROOT.parents[2] / "mevo_fast_pzt_f10o_20261013/computation_examples/mevo_acoustic_microbubbles"
    assert (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=original, text=True).strip()
        == p["baseline_commit"]
    )
    assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=original, text=True).strip()
    for rel, r in read(original / "results/straight_tube_exp1/f10o/ARTIFACT_MANIFEST.json")["files"].items():
        assert sha(original / rel) == r["sha256"], rel
    for rel in [
        "fast_sim/f10g_contact.py",
        "fast_sim/f10h_contact.py",
        "fast_sim/water_b2.py",
        "fast_sim/water_control_engine.py",
        "fast_sim/pdms_cavity_engine.py",
        "fast_sim/cubic_field.py",
        "acoustics/water_disc_proxy.py",
        "coupling/pdms_cavity.py",
    ]:
        assert sha(ROOT / "src/mevo_mb" / rel) == sha(original / "src/mevo_mb" / rel), rel
    checks["F10O_original_clean_all_historical_artifacts_and_protected_models_unchanged"] = True
    reg = read(BASE / "PREREGISTRATION.json")
    assert reg["protocol"] == p
    for rel, h in reg["input_sha256"].items():
        assert sha(OLD / rel) == h, rel
    assert p["XY_targets_m"] == {"X": [2e-7, 0, 0], "Y": [0, 2e-7, 0]}
    subprocess.run(["git", "merge-base", "--is-ancestor", "546f286f1", "36343716a"], cwd=ROOT, check=True)
    subprocess.run(["git", "merge-base", "--is-ancestor", "36343716a", "cc35fad01"], cwd=ROOT, check=True)
    checks["cleanup_and_thresholds_committed_before_science_exact_inputs_preserved"] = True
    for rel, r in read(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert sha(BASE / rel) == r["sha256"], rel
    checks["all_remote_artifacts_individually_verified"] = True
    initial = np.load(OLD / "initial_states.npz")["aggregate_m"]
    count = 0
    failures = []
    for f in BASE.glob("*/CAMPAIGN.json"):
        c = read(f)
        ident = c["identity"]
        dep = read(BASE / "execution" / (f.parent.name + ".json"))
        assert ident["execution_commit"] == dep["sha"]
        pack = BASE / "execution" / (dep["sha"] + ".tar.gz")
        assert sha(pack) == dep["archive_sha256"]
        with tarfile.open(pack) as arc:
            for rel, h in ident["source_sha256"].items():
                assert hashlib.sha256(arc.extractfile(rel).read()).hexdigest() == h, rel
        for result in (f.parent / "runs").glob("*/RESULT.json"):
            count += 1
            r = read(result)
            folder = result.parent
            inp = read(folder / "INPUT.json")
            assert inp["identity"]["execution_commit"] == ident["execution_commit"]
            assert inp["preregistration_sha256"] == sha(BASE / "PREREGISTRATION.json")
            old = read(OLD / "campaign_v1/runs" / inp["source_F10N"] / "INPUT.json")["events"][0]
            for k in ["amplitude_Pa", "phase_rad", "enabled", "flow_enabled"]:
                assert inp["control"][k] == old[k]
            assert sha(folder / "trajectory.h5") == r["trajectory_sha256"]
            with h5py.File(folder / "trajectory.h5") as h:
                t = h["time_s"][:]
                x = h["positions_m"][:]
                m = h["step_metrics"][:]
                accepted = m[:, 0].astype(bool)
                np.testing.assert_array_equal(x[0], initial)
                np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
                assert np.isfinite(x).all() and h.attrs["flow"] == "STRICT_ZERO"
                assert len(x) == r["accepted_steps"] + 1 and len(m) == r["attempted_steps"]
                assert np.all(abs(np.diff(t) - inp["dt_s"]) < 1e-14)
                # All these center-region states lie strictly inside the finite main lumen;
                # this independently proves no PDMS/closed-X penetration for these data.
                assert np.all(abs(x[:, :, 0]) + p["radius_m"] < 0.01625)
                assert np.all(np.linalg.norm(x[:, :, 1:], axis=-1) + p["radius_m"] < 0.0008)
                actual = np.array([max(0, -np.min(pdist(row) - 2 * p["radius_m"])) for row in x])
                assert actual.max() <= p["hard_geometry_m"]
                for col, lim in [
                    (1, p["mobility_residual_m"]),
                    (2, p["hard_geometry_m"]),
                    (4, p["contact_correction_m"]),
                    (12, p["primal_violation_m"]),
                    (13, p["dual_negative_N"]),
                    (15, p["complementarity_equivalent_gap_m"]),
                    (17, p["actual_geometry_complementarity_m"]),
                ]:
                    assert np.all(m[accepted, col] <= lim), (folder.name, col)
                assert np.count_nonzero(m[:, 5]) == 0
                pressure = abs(h["force_samples/pressure_Pa"][:])
                assert np.all((pressure == 0) | ((pressure >= 0.1) & (pressure <= 40000)))
                if r["status"] == "PASS":
                    assert accepted.all() and abs(t[-1] - inp["duration_s"]) < 1e-13
                else:
                    failures.append((f.parent.name, folder.name, r["reasons"]))
                    assert not accepted[-1] and accepted[:-1].all() and "failure_proposal_m" in h
                    assert not np.array_equal(x[-1], h["failure_proposal_m"][:])
                if r["source"] == "OFF":
                    np.testing.assert_array_equal(x, np.broadcast_to(x[0], x.shape))
    checks["immutable_executions_true_geometry_budgets_fixed_ids_zero_flow_and_failure_prefixes"] = True
    assert len(failures) == 3 and count == 12, (count, failures)
    costs = read(BASE / "analysis/COSTS.json")
    assert costs["total_numerical_s"] < 600
    g = read(BASE / "analysis/FINAL_GATE.json")
    assert not g["navigation_admitted"] and not g["training_ready"]
    assert g["groups"]["STRONG"]["BASIC_NUMERICS_PASS"] == "FAIL"
    assert (
        g["groups"]["WEAK"]["simulation_group_position_eligible"]
        and g["groups"]["WEAK"]["simulation_group_shape_eligible"]
    )
    # A synthetic metadata fixture tests gate independence, not a scientific success case.
    rr = read(BASE / "analysis/RUN_ANALYSIS.json")["WEAK_DT1"]
    cc = read(BASE / "analysis/COMPARISONS.json")
    d = copy.deepcopy(cc["WEAK_DT05"])
    ca = copy.deepcopy(cc["WEAK_FINE"])
    d["fixed_id_m"] = 1e-6
    test = layered(rr, d, ca, [rr["topology"]] * 3, p)
    assert (
        test["simulation_group_position_eligible"]
        and test["simulation_group_shape_eligible"]
        and not test["simulation_particle_id_eligible"]
    )
    checks["group_ID_gates_independent_and_no_invalid_navigation_promotion"] = True
    for version in ["v1", "v2"]:
        assert read(BASE / ("local_foundations_" + version + ".json"))["status"] == "PASS"
    checks["small_independent_references_and_fail_closed_fixtures"] = True
    if final:
        prov = read(BASE / "figures/FIGURE_PROVENANCE.json")
        review = read(BASE / "figures/VISUAL_REVIEW.json")
        assert len(prov["figures"]) == 4 and review["status"] == "PASS"
        assert sha(ROOT / "scripts/straight_tube_exp1/f10p/figures.py") == prov["source_script_sha256"]
        for rel, h in prov["input_sha256"].items():
            assert sha(ROOT / rel) == h, rel
        for rec in review["figures"]:
            assert rec["opened_png"] and rec["opened_pdf_page_preview"]
            for ext, h in rec["file_sha256"].items():
                assert sha(BASE / "figures" / (rec["figure"] + "." + ext)) == h
        checks["four_actual_PNG_PDF_figures_and_visual_review"] = True
        d = BASE / "progress_document"
        before = (d / "before.md").read_bytes()
        after = (d / "after.md").read_bytes()
        ch = read(d / "CHECK.json")
        assert sha(ch["path"]) == ch["after_sha256"] and sha(ch["backup"]) == ch["before_sha256"]
        i = after.index("\n**F10P ".encode())
        j = after.index("\n---\n\n### 尾声".encode(), i)
        assert after[:i] + after[j:] == before
        assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
        checks["progress_byte_preservation_backup_hashes_and_heading_structure"] = True
    return dict(
        status="PASS",
        checks=checks,
        run_count=count,
        failures=failures,
        scope="ENGINEERING evidence integrity only; strong full-window BASIC remains FAIL",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", action="store_true")
    args = parser.parse_args()
    r = verify(args.final)
    write(
        BASE / "analysis" / ("FINAL_VERIFICATION.json" if args.final else "ENGINEERING_VERIFICATION.json"), r
    )
    print(r)


if __name__ == "__main__":
    main()

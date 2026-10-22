"""Read-only evidence verification. --record additionally writes local manifests."""

import argparse
import re
from pathlib import Path
import numpy as np
import h5py
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .common import BASE, ROOT, CONFIG, configs, read, write, sha


def manifest():
    files = [
        p
        for p in BASE.rglob("*")
        if p.is_file() and p.name not in ["ARTIFACT_MANIFEST.json", "VERIFICATION.json"]
    ]
    files += list(CONFIG.glob("*.yaml")) + list((ROOT / "scripts/straight_tube_exp1/f10n").glob("*.py"))
    files += [
        ROOT / p
        for p in [
            "src/mevo_mb/acoustics/water_disc_proxy.py",
            "src/mevo_mb/fast_sim/cubic_field.py",
            "src/mevo_mb/fast_sim/water_control_engine.py",
            "docs/F10N_CONTROL_MODEL.md",
            "docs/FAST_PZT_F10N_DEVELOPMENT.md",
        ]
    ]
    data = {str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(set(files))}
    write(BASE / "ARTIFACT_MANIFEST.json", dict(files=data, bytes=sum(r["bytes"] for r in data.values())))


def verify():
    checks = {}
    c = configs()
    p = c["protocol"]
    r = p["radius_m"]
    reg = read(BASE / "PREREGISTRATION.json")
    assert read(BASE / "CONFIG_SNAPSHOT.json") == c and reg["protocol"] == p
    assert reg["target_locked_before_any_GPU_trajectory"] and p["target_position_m"] == [3e-6, 2e-6, 1e-6]
    assert sha(BASE / "initial_states.npz") == reg["initial_sha256"]
    for key in ["coarse", "fine"]:
        assert sha(BASE / "cache" / f"{key}.h5") == reg["cache_sha256"][key]
    assert p["flow"] == "STRICT_ZERO" and not p["training_ready"]
    checks["registered_inputs_targets_and_cache_hashes_unchanged"] = True
    cpu = read(BASE / "analysis/CPU_CHECKS.json")
    assert cpu["status"] == "PASS" and all(cpu["checks"].values())
    accepted = read(BASE / "ACCEPTED_CAMPAIGN.json")
    out = BASE / accepted["path"]
    assert sha(out / "CAMPAIGN.json") == accepted["sha256"]
    campaign = read(out / "CAMPAIGN.json")
    gpu = read(out / "GPU_CHECKS.json")
    assert gpu["status"] == "PASS" and all(gpu["checks"].values())
    for rel, digest in campaign["identity"]["files_sha256"].items():
        if rel.startswith("src/") or rel.endswith(("campaign.py", "preflight.py")):
            assert sha(ROOT / rel) == digest, rel
    checks["CPU_integral_and_actual_GPU_kernel_checks"] = True
    # Four independently gated sources, checked against actual saved basis fields.
    basis = CachedWaterBasis(BASE / "cache/coarse.h5")
    x = np.load(BASE / "initial_states.npz")["aggregate_m"]
    pb, gb = basis.query_basis(x)
    for j in range(4):
        enabled = np.zeros(4, bool)
        enabled[j] = True
        u = basis.vector([100] * 4, [0.2, 0.4, 0.6, 0.8], enabled)
        assert np.count_nonzero(u) == 1
        np.testing.assert_allclose(u @ pb, 100 * np.exp(1j * [0.2, 0.4, 0.6, 0.8][j]) * pb[j], atol=1e-10)
    np.testing.assert_array_equal(basis.vector([100] * 4, [0] * 4, [False] * 4), np.zeros(4))
    checks["independent_source_enable_and_phase"] = True
    trials = read(out / "TRIALS.json")
    assert len(trials) == campaign["trajectory_count"] <= p["max_total_trajectories"]
    assert campaign["wall_s"] < p["campaign_wall_budget_s"]
    failed = []
    for row in trials:
        folder = out / "runs" / row["name"]
        inp = read(folder / "INPUT.json")
        summary = read(folder / "SUMMARY.json")
        assert inp["preregistration_sha256"] == sha(BASE / "PREREGISTRATION.json")
        assert sha(folder / "trajectory.h5") == row["trajectory_sha256"]
        assert inp["events"][0]["duration_s"] == p["duration_s"] and not inp["events"][0]["flow_enabled"]
        with h5py.File(folder / "trajectory.h5") as h:
            x = h["positions_m"][:, 0]
            t = h["time_s"][:]
            pp = abs(h["pressure_Pa"][:])
            np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
            assert np.isfinite(x).all() and np.all(np.diff(t) > 0)
            assert np.all(abs(x[:, :, 0]) + r <= 0.01625) and np.all(
                np.linalg.norm(x[:, :, 1:], axis=-1) + r < 0.0008
            )
            assert np.all((pp == 0) | ((pp >= 0.1) & (pp <= 40000)))
            for key in ["flow/velocity_m_s", "flow/gradient_s_inv"]:
                assert key in h and np.count_nonzero(h[key][:]) == 0
            if row["name"] in ["OFF", "DISPERSED_OFF"]:
                np.testing.assert_array_equal(x, np.broadcast_to(x[0], x.shape))
                for key in ["pressure_Pa", "primary_force_N", "secondary_force_N"]:
                    assert np.count_nonzero(h[key][:]) == 0
            if summary["status"] == "FAILED":
                failed.append(row["name"])
                assert "failure_proposal" in h and summary["reason"]
            else:
                assert t[-1] == p["duration_s"] and row["eligible"]
                q = row["quality"]
                assert q["error_bits"] == 0 and q["max_accepted_penetration_m"] <= 0.01 * r
                assert (
                    q["max_contact_pass_correction_m"] <= 0.1 * r
                    and q["maximum_wall_correction_m"] <= 0.1 * r
                )
            np.testing.assert_allclose(
                x[-1].mean(0) - x[0].mean(0), row["displacement_m"], rtol=0, atol=1e-18
            )
    assert failed == ["NAV_STATIC_0"]
    checks["all_21_actual_trajectories_fixed_ID_flow_pressure_walls_and_failure"] = True
    gate = read(BASE / "analysis/FINAL_GATE.json")
    sens = read(BASE / "analysis/NUMERICAL_SENSITIVITY.json")
    assert not gate["target_reached"] and not gate["shape_target_reached"]
    assert (
        not gate["navigation_particle_resolution_pass"]
        and gate["navigation_fixed_id_numerical_envelope_um"] > 6
    )
    assert sens["navigation_dt"]["max_fixed_id_position_difference_m"] > p["position_dt_tolerance_m"]
    assert gate["shape_directional_effect_exceeds_tested_numerical_difference"] and not gate["training_ready"]
    checks["no_target_success_or_particle_convergence_claim"] = True
    remote = read(BASE / "REMOTE_MANIFEST.json")
    for rel, rec in remote["files"].items():
        assert sha(BASE / rel) == rec["sha256"], rel
    checks["all_remote_files_returned_with_exact_hashes"] = True
    prov = read(BASE / "figures/FIGURE_PROVENANCE.json")
    review = read(BASE / "figures/VISUAL_REVIEW.json")
    assert len(prov["figures"]) == len(review["figures"]) == 4 and review["status"] == "PASS"
    assert sha(ROOT / "scripts/straight_tube_exp1/f10n/figures.py") == prov["source_script_sha256"]
    for rec in prov["figures"]:
        for ext, digest in rec["file_sha256"].items():
            assert sha(BASE / "figures" / (rec["figure"] + "." + ext)) == digest
        for rel, digest in rec["source_sha256"].items():
            assert sha(ROOT / rel) == digest, rel
    for rec in review["figures"]:
        assert rec["opened_png"] and rec["opened_pdf_page_preview"]
        for ext, digest in rec["file_sha256"].items():
            assert sha(BASE / "figures" / (rec["figure"] + "." + ext)) == digest
    checks["four_actual_reviewed_PNG_PDF_figures"] = True
    pg = read(BASE / "progress_document/CHECK.json")
    before = (BASE / "progress_document/before.md").read_bytes()
    after = (BASE / "progress_document/after.md").read_bytes()
    assert sha(pg["path"]) == pg["after_sha256"] and sha(pg["backup"]) == pg["before_sha256"]
    i = after.index("\n**F10N ".encode())
    j = after.index("\n---\n\n### 尾声".encode(), i)
    assert after[:i] + after[j:] == before
    assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
    checks["progress_backup_one_paragraph_original_bytes_preserved"] = True
    history = {}
    for stage in ["f10l", "f10m", "f10m_r3"]:
        root = (
            Path("/home/lzy/projects/computation_examples/_archives")
            / ("mevo_fast_pzt_" + stage + "_20261013")
            / "computation_examples/mevo_acoustic_microbubbles"
        )
        path = root / "results/straight_tube_exp1" / stage / "ARTIFACT_MANIFEST.json"
        items = read(path)["files"]
        for rel, rec in items.items():
            assert sha(root / rel) == rec["sha256"], str(root / rel)
        history[stage] = dict(status="PASS", unchanged_files=len(items))
    if (BASE / "ARTIFACT_MANIFEST.json").exists():
        for rel, rec in read(BASE / "ARTIFACT_MANIFEST.json")["files"].items():
            assert sha(ROOT / rel) == rec["sha256"], rel
        checks["local_artifact_manifest"] = True
    return dict(
        status="PASS",
        checks=checks,
        historical_workspaces=history,
        scientific_status="Both targets unmet; navigation particle morphology unresolved; uncalibrated weak shape response measured",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--record", action="store_true")
    a = p.parse_args()
    if a.record:
        manifest()
    result = verify()
    if a.record:
        write(BASE / "VERIFICATION.json", result)
    print(result)

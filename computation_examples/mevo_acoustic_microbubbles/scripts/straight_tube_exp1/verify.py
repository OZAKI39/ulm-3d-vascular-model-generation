"""Read-only full F10L delivery verification; no solver reruns."""

import argparse
import hashlib
import json
import re
import h5py
import numpy as np
from mevo_mb.coupling.straight_tube import apparatus
from .common import BASE, ROOT, read, sha, write, validate_saved_acoustic_inputs


def run():
    checks = {}
    assert validate_saved_acoustic_inputs()
    config = read(BASE / "CONFIG_SNAPSHOT.json")
    assert config["lab_parameters"]["parameters"]["wall_thickness_m"]["value"] is None
    assert config["demo"]["flow"] == {"mode": "zero"}
    tube, elements = apparatus(config["demo"]["geometry"], config["acoustic"]["acoustics"]["frequency_Hz"])
    assert len(elements) == 3 and all(e["size_m"] == [0.009] for e in elements)
    pre = read(BASE / "analysis/PREFLIGHT.json")
    assert pre["geometry"]["status"] == "PASS" and pre["geometry"]["no_local_z_sources"]
    assert all(e["gap_matches_registered"] for e in pre["geometry"]["pzt"])
    checks["three_9mm_faces_provisional_ID_actual_gaps_and_unknown_lab_wall"] = True
    response = ROOT / config["bubbles"]["response_path"]
    assert sha(response) == config["bubbles"]["response_sha256"]
    checks["original_response_identity_and_no_historical_field_reuse"] = True
    source = read(BASE / "acoustic/SOURCE_GRID.json")
    assert source["source_count"] == 3 and min(source["edge_clearance_cells"]) >= 1 and source["pml_outside"]
    assert source["source_weights_sha256"] == sha(BASE / "acoustic/source_weights.npz")
    basis = BASE / "acoustic/basis.h5"
    ac = read(BASE / "acoustic/ACOUSTIC_GATE.json")
    assert ac["status"] == "PASS" and ac["basis_sha256"] == sha(basis)
    assert ac["whole_acoustic_s"] < 600 and ac["finite_difference_gradient_relative_error"] < 1e-6
    with h5py.File(basis) as h:
        p = h["P_complex"][:]
        assert p.shape[0] == 3 and np.isfinite(p).all() and h["valid_roi_mask"][:].all()
        meta = json.loads(h.attrs["metadata_json"])
        assert "F6 exact derivative" in meta["gradient"]
        assert h.attrs["harmonic_convention"] == "Re[P exp(+i omega t)]"
        for i in range(3):
            m = read(BASE / f"acoustic/basis_{i:02d}/metadata.json")
            assert m["recording_half_window_relative_L2"] < 0.01
            assert m["sampled_device_peak_VRAM_MiB"] < 16 * 1024
            control = np.array(m["source_control_real_Pa"])
            assert np.array_equal(control, np.eye(3)[i] * 1000)
    checks["three_actual_complex_bases_source_support_phasor_gradient_and_budget"] = True
    engine = read(BASE / "analysis/ENGINE_CHECKS.json")
    assert engine["status"] == "PASS" and all(engine["checks"].values())
    assert len(engine["checks"]) >= 21
    checks["21_CPU_GPU_guards_and_bitwise_M0_equivalence"] = True
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    R = config["bubbles"]["radius_m"]
    dg = read(BASE / "dynamics/DYNAMICS_GATE.json")
    for name in ("OFF", "THREE_PZT_ON"):
        folder = BASE / "dynamics" / name
        tr = np.load(folder / "trajectory.npz")
        x = tr["positions_m"]
        assert x.shape == (101, 1, 32, 3) and tube.contains(x, R).all()
        np.testing.assert_array_equal(tr["particle_id"], np.arange(32))
        np.testing.assert_allclose(tr["time_s"], np.arange(101) * 1e-5, rtol=0, atol=2e-18)
        obs = np.load(folder / "frame_observables.npz")
        assert np.count_nonzero(obs["background_velocity_m_s"]) == 0
        assert np.isfinite(obs["pressure_complex_Pa"]).all()
        if name == "OFF":
            np.testing.assert_array_equal(x, np.broadcast_to(initial, x.shape))
            assert np.count_nonzero(obs["primary_N"]) == np.count_nonzero(obs["secondary_N"]) == 0
        assert dg["cases"][name]["quality"]["error_bits"] == 0
        assert dg["cases"][name]["quality"]["max_accepted_penetration_m"] <= 0.01 * R
        assert dg["cases"][name]["quality"]["max_contact_pass_correction_m"] <= 0.1 * R
        assert dg["cases"][name]["quality"]["maximum_wall_proposal_m"] <= 0.1 * R
    assert read(BASE / "dynamics/NO_CFD_ACCESS.json")["CFD_access_count"] == 0
    checks["two_real_1ms_32ID_trajectories_OFF_exact_no_drift_zero_flow_no_CFD"] = True
    rem = read(BASE / "REMOTE_MANIFEST.json")
    for rel, r in rem["files"].items():
        assert sha(BASE / rel) == r["sha256"] and (BASE / rel).stat().st_size == r["bytes"]
    checks["all_remote_scientific_files_and_failure_logs_hash_verified"] = True
    provenance = read(BASE / "figures/FIGURE_PROVENANCE.json")
    review = read(BASE / "figures/VISUAL_REVIEW.json")
    assert len(provenance["figures"]) == len(review["figures"]) == 4 and review["status"] == "PASS"
    for fig in provenance["figures"]:
        for ext, digest in fig["file_sha256"].items():
            assert sha(BASE / "figures" / (fig["figure"] + "." + ext)) == digest
        for rel, digest in fig["source_sha256"].items():
            assert sha(ROOT / rel) == digest
    assert sha(ROOT / "scripts/straight_tube_exp1/figures.py") == provenance["source_script_sha256"]
    checks["four_PNG_PDF_pairs_actual_visual_review_and_provenance"] = True
    prog = read(BASE / "progress_document/CHECK.json")
    before = (BASE / "progress_document/before.md").read_bytes()
    after = (BASE / "progress_document/after.md").read_bytes()
    assert sha(prog["backup"]) == prog["before_sha256"] and sha(prog["path"]) == prog["after_sha256"]
    assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
    marker = b"\n**F10L "
    assert after.count(marker) == 1
    start = after.index(marker)
    end = after.index("\n---\n\n### 尾声".encode(), start)
    assert after[:start] + after[end:] == before
    checks["one_progress_paragraph_backup_all_old_bytes_headings_preserved"] = True
    gate = read(BASE / "analysis/FINAL_GATE.json")
    assert (
        gate["engineering_status"] == "PASS" and not gate["training_ready"] and not gate["experiment_matched"]
    )
    assert set(gate["statuses"]) == {
        "EXPERIMENTAL_GEOMETRY_PENDING",
        "MATERIAL_ACOUSTICS_UNCALIBRATED",
        "SOURCE_PRESSURE_UNCALIBRATED",
        "PHYSICS_UNCALIBRATED",
    }
    checks["uncalibrated_science_no_navigation_or_AI_admission"] = True
    manifest = BASE / "ARTIFACT_MANIFEST.json"
    if manifest.exists():
        for rel, r in read(manifest)["files"].items():
            assert sha(ROOT / rel) == r["sha256"]
        checks["local_complete_artifact_manifest"] = True
    return dict(status="PASS", checks=checks)


def manifest():
    paths = [
        p
        for p in BASE.rglob("*")
        if p.is_file() and p.name not in ("ARTIFACT_MANIFEST.json", "VERIFICATION.json")
    ]
    paths += list((ROOT / "configs/straight_tube_exp1").glob("*.yaml"))
    paths += list((ROOT / "scripts/straight_tube_exp1").glob("*.py"))
    paths += [
        ROOT / "src/mevo_mb/coupling/straight_tube.py",
        ROOT / "src/mevo_mb/fast_sim/straight_tube_engine.py",
        ROOT / "src/mevo_mb/fast_sim/f3_engine.py",
    ]
    paths += [
        ROOT / "docs" / name
        for name in (
            "FAST_PZT_F10L_STRAIGHT_TUBE_DEVELOPMENT.md",
            "STRAIGHT_TUBE_EXPERIMENT_MODEL.md",
            "STRAIGHT_TUBE_LAB_PARAMETERS.md",
        )
    ]
    files = {
        str(p.relative_to(ROOT)): dict(
            bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()
        )
        for p in sorted(paths)
    }
    write(BASE / "ARTIFACT_MANIFEST.json", dict(files=files, bytes=sum(r["bytes"] for r in files.values())))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    if args.record:
        manifest()
    report = run()
    if args.record:
        write(BASE / "VERIFICATION.json", report)
    print(json.dumps(report, indent=2))

"""Read-only delivery verification; --record also stores the local manifest."""

import argparse
import ast
import re
from pathlib import Path
import numpy as np
import h5py
from mevo_mb.coupling.pdms_cavity import apparatus
from .common import BASE, ROOT, CONFIG, configs, read, write, sha


def verify():
    checks = {}
    c = configs()
    assert read(BASE / "CONFIG_SNAPSHOT.json") == c
    cavity, elements = apparatus(c["geometry"])
    assert len(elements) == 4 and cavity.lumen.ends == ("closed", "closed")
    assert np.array_equal(cavity.ports, [-0.01085, 0.01085])
    assert cavity.a == 0.00075 and cavity.R == 0.0008
    assert all(
        p["reference_status"] == "PROVISIONAL_MEASUREMENT_REFERENCE" for p in c["geometry"]["top_ports"]
    )
    checks["confirmed_four_faces_two_open_tops_two_sealed_X"] = True
    assert read(BASE / "analysis/CONNECTIVITY.json")["components"] == 1
    cpu = read(BASE / "analysis/CPU_CHECKS.json")
    assert cpu["status"] == "PASS" and all(cpu["checks"].values())
    checks["CPU_geometry_connectivity_and_finite_radius_checks"] = True
    material = BASE / "materials/regions_and_unresolved_media.h5"
    assert sha(material) == read(BASE / "materials/IDENTITY.json")["sha256"]
    with h5py.File(material) as h:
        axes = [h[k + "_m"][:] for k in "xyz"]
        # Independent evaluation of sampled saved labels against current geometry.
        points = np.stack(np.meshgrid(axes[0][::23], axes[1][::11], axes[2][::11], indexing="ij"), -1)
        np.testing.assert_array_equal(h["region_id"][::23, ::11, ::11], cavity.labels(points))
        np.testing.assert_array_equal(h["boundary_id"][::23, ::11, ::11], cavity.boundary_labels(points))
        iy = np.argmin(abs(axes[1]))
        iz = np.argmin(abs(axes[2] - 0.003))
        for i, x in enumerate(cavity.ports):
            ix = np.argmin(abs(axes[0] - x))
            assert h["region_id"][ix, iy, iz] == 3 + i and np.isnan(h["sound_speed_map"][ix, iy, iz])
        for i, x in enumerate([-cavity.L, cavity.L]):
            ix = np.argmin(abs(axes[0] - x))
            iz0 = np.argmin(abs(axes[2]))
            assert h["boundary_id"][ix, iy, iz0] == i + 1 and np.isnan(h["sound_speed_map"][ix, iy, iz0])
    checks["unknown_phase_and_seal_acoustics_not_filled_with_water"] = True
    accepted = read(BASE / "ACCEPTED_VALIDATION.json")
    out = BASE / accepted["path"]
    assert sha(out / "VALIDATION.json") == accepted["validation_sha256"]
    val = read(out / "VALIDATION.json")
    assert val["status"] == "PASS" and all(val["checks"].values())
    for rel in [
        "src/mevo_mb/coupling/pdms_cavity.py",
        "src/mevo_mb/fast_sim/pdms_cavity_engine.py",
        "scripts/straight_tube_exp1/f10m_r3/validate_gpu.py",
    ]:
        assert sha(ROOT / rel) == val["identity"]["files_sha256"][rel]
    wall = read(out / "WALL_CHECKS.json")
    assert wall["distance_max_error_m"] < c["protocol"]["wall_query_tolerance_m"]
    assert wall["projection_max_error_m"] < 5e-10
    source = read(out / "source_input/CHECK.json")
    assert source["status"] == "PASS" and source["independent_weight_rank"] == 4
    assert source["actual_wave_solves"] == 0
    checks["current_GPU_wall_code_and_four_source_wiring_validated"] = True
    tree = ast.parse((ROOT / "src/mevo_mb/fast_sim/pdms_cavity_engine.py").read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "PDMSCavitySimulator")
    assert [n.id for n in cls.bases] == ["FullParticleBatch"]
    imported = [a.name for n in tree.body if isinstance(n, ast.ImportFrom) for a in n.names]
    assert "StraightTubeSimulator" not in imported and "cylinder_project" not in imported
    checks["R3_does_not_inherit_legacy_open_endpoint_adapter"] = True
    with h5py.File(out / "OFF_trajectory.h5") as h:
        x = h["positions_m"][:]
        assert x.shape == (11, 1, 16, 3) and np.isclose(h["time_s"][-1], 0.0001)
        np.testing.assert_array_equal(h["particle_id"][:], np.arange(16))
        np.testing.assert_array_equal(x, np.broadcast_to(x[0], x.shape))
        assert all(cavity.query(q, c["protocol"]["bubbles"]["radius_m"])["valid"] for q in x[0, 0])
        assert np.count_nonzero(h["primary_force_N"][:]) == np.count_nonzero(h["secondary_force_N"][:]) == 0
        assert h.attrs["end_conditions"] == "closed,closed"
    off = read(out / "OFF_SUMMARY.json")
    assert off["status"] == "COMPLETE" and off["quality"]["error_bits"] == 0
    assert off["quality"]["maximum_wall_correction_m"] == off["quality"]["max_accepted_penetration_m"] == 0
    assert read(out / "NO_CFD_ACCESS.json")["CFD_access_count"] == 0
    checks["actual_16_ID_0p1ms_OFF_static_legal_no_CFD"] = True
    manifest = read(BASE / "REMOTE_MANIFEST.json")
    for rel, rec in manifest["files"].items():
        p = BASE / rel
        assert p.stat().st_size == rec["bytes"] and sha(p) == rec["sha256"], rel
    assert (BASE / "validation_v1/FAILURE.json").exists()
    assert (BASE / "failed_attempts/preflight_v1/stderr.log").exists()
    checks["all_remote_hashes_and_failed_attempts_retained"] = True
    prov = read(BASE / "figures/FIGURE_PROVENANCE.json")
    review = read(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and len(review["figures"]) == len(prov["figures"]) == 4
    assert sha(ROOT / "scripts/straight_tube_exp1/f10m_r3/figures.py") == prov["source_script_sha256"]
    for f in prov["figures"]:
        for ext, digest in f["file_sha256"].items():
            assert sha(BASE / "figures" / (f["figure"] + "." + ext)) == digest
        for rel, digest in f["source_sha256"].items():
            assert sha(ROOT / rel) == digest
    for f in review["figures"]:
        for ext, digest in f["file_sha256"].items():
            assert sha(BASE / "figures" / (f["figure"] + "." + ext)) == digest
        assert f["opened_png"] and f["opened_pdf_page_preview"]
    checks["four_reviewed_actual_PNG_PDF_pairs_and_sources"] = True
    pg = read(BASE / "progress_document/CHECK.json")
    before = (BASE / "progress_document/before.md").read_bytes()
    after = (BASE / "progress_document/after.md").read_bytes()
    assert sha(pg["path"]) == pg["after_sha256"] and sha(pg["backup"]) == pg["before_sha256"]
    a = after.index("\n**F10M-R3 ".encode())
    b = after.index("\n---\n\n### 尾声".encode(), a)
    assert after[:a] + after[b:] == before
    assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
    checks["progress_backup_one_paragraph_all_old_bytes_preserved"] = True
    gate = read(BASE / "analysis/FINAL_GATE.json")
    assert (
        gate["engineering_status"] == "PASS" and not gate["training_ready"] and not gate["experiment_matched"]
    )
    assert gate["full_apparatus_fields"] == gate["ON"] == "NOT_RUN" and gate["actual_wave_solves"] == 0
    assert not list(BASE.rglob("field.npz"))
    checks["no_invented_apparatus_field_ON_or_surface_calibration"] = True
    if (BASE / "ARTIFACT_MANIFEST.json").exists():
        for rel, rec in read(BASE / "ARTIFACT_MANIFEST.json")["files"].items():
            assert sha(ROOT / rel) == rec["sha256"], rel
        checks["local_complete_manifest"] = True
    # Original workspaces, when available, remain independently auditable.
    historical = {}
    for stage in ["f10l", "f10m"]:
        root = (
            Path("/home/lzy/projects/computation_examples/_archives")
            / ("mevo_fast_pzt_" + stage + "_20261013")
            / "computation_examples/mevo_acoustic_microbubbles"
        )
        record = root / ("results/straight_tube_exp1/" + stage + "/ARTIFACT_MANIFEST.json")
        if record.exists():
            items = read(record)["files"]
            for rel, rec in items.items():
                assert sha(root / rel) == rec["sha256"], str(root / rel)
            historical[stage] = dict(status="PASS", unchanged_files=len(items))
    return dict(status="PASS", checks=checks, historical_workspaces=historical)


def manifest():
    files = [
        p
        for p in BASE.rglob("*")
        if p.is_file() and p.name not in ["ARTIFACT_MANIFEST.json", "VERIFICATION.json"]
    ]
    files += list(CONFIG.glob("*.yaml")) + list((ROOT / "scripts/straight_tube_exp1/f10m_r3").glob("*.py"))
    files += [
        ROOT / p
        for p in [
            "src/mevo_mb/coupling/pdms_cavity.py",
            "src/mevo_mb/fast_sim/pdms_cavity_engine.py",
            "docs/FAST_PZT_F10M_R3_DEVELOPMENT.md",
            "docs/PDMS_FOUR_PZT_OPEN_TOP_SEALED_ENDS.md",
        ]
    ]
    data = {str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(set(files))}
    write(BASE / "ARTIFACT_MANIFEST.json", dict(files=data, bytes=sum(r["bytes"] for r in data.values())))


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

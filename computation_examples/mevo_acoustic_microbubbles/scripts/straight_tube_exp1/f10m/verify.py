"""Read-only final artifact checks; --record writes the complete SHA manifest."""

import argparse
import re
import numpy as np
import h5py
from .common import BASE, ROOT, CONFIG, read, write, sha, configs
from mevo_mb.coupling.pdms_block import apparatus, material_labels


def run():
    checks = {}
    c = configs()
    assert read(BASE / "CONFIG_SNAPSHOT.json") == c
    g = read(BASE / "analysis/GEOMETRY.json")
    lumen, src = apparatus(c["geometry"])
    assert (
        g["block_dimensions_m"] == [0.0325, 0.009, 0.009]
        and g["lumen_radius_m"] == 0.0008
        and g["centered"]
        and not g["independent_tube_wall"]
    )
    assert len(src) == 3 and all(e["size_m"] == [0.007] for e in src)
    assert np.isclose(g["side_PDMS_path_m"], 0.0037) and g["end_opening_coverage_fraction"] == 1
    assert (
        g["top_port"]["diameter_m"] == 0.0015
        and g["top_port"]["boundary"] == "open"
        and g["end_conditions"] == ["open_fail_stop", "closed"]
    )
    checks["confirmed_PDMS_7mm_sources_top_port_and_closed_X_plus"] = True
    m = BASE / "materials/material_maps.h5"
    assert sha(m) == read(BASE / "materials/IDENTITY.json")["sha256"]
    with h5py.File(m) as h:
        assert h["sound_speed_map"].shape == tuple(read(BASE / "analysis/RESOURCE_PREFLIGHT.json")["shape"])
        ix = np.argmin(abs(h["x_m"][:] - 0.01085))
        iy = np.argmin(abs(h["y_m"][:]))
        iz = np.argmin(abs(h["z_m"][:] - 0.003))
        assert h["material_id"][ix, iy, iz] == 2
        coords = np.stack(np.meshgrid(h["x_m"][::67], h["y_m"][::23], h["z_m"][::23], indexing="ij"), axis=-1)
        np.testing.assert_array_equal(
            h["material_id"][::67, ::23, ::23], material_labels(coords, c["geometry"])
        )
    checks["hashed_three_material_maps_include_top_port"] = True
    for rel in ["analysis/CPU_CHECKS.json", "adapter_validation/VALIDATION.json", "source_bridge/GATE.json"]:
        d = read(BASE / rel)
        assert d["status"] == "PASS"
        if "checks" in d:
            assert all(d["checks"].values())
    assert read(BASE / "adapter_validation/VALIDATION.json")["frequency_Hz"] == 5e6
    assert read(BASE / "source_bridge/GATE.json")["scope"].startswith("MINIATURE")
    checks["bounded_GPU_heterogeneous_validation_separate_from_apparatus"] = True
    assert not (BASE / "acoustic/basis.h5").exists()
    assert read(BASE / "analysis/ACOUSTIC_ADMISSION.json")["status"] == "NOT_RUN"
    assert read(BASE / "dynamics/ON_NOT_RUN.json")["status"] == "NOT_RUN"
    checks["no_fabricated_apparatus_field_or_ON"] = True
    d = read(BASE / "dynamics/SUMMARY.json")
    assert d["status"] == "COMPLETE" and all(d["checks"].values())
    assert d["end_conditions"] == ["open_fail_stop", "closed"]
    with h5py.File(BASE / "dynamics/OFF_trajectory.h5") as h:
        x = h["positions_m"][:]
        np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
        assert x.shape == (101, 1, 32, 3) and lumen.contains(x, 1.78e-6).all()
        np.testing.assert_array_equal(x, np.broadcast_to(x[0], x.shape))
        for key in ["pressure_complex_Pa", "primary_force_N", "secondary_force_N"]:
            assert np.count_nonzero(h[key][:]) == 0
        assert np.isclose(h["time_s"][-1], 0.001)
    assert read(BASE / "dynamics/NO_CFD_ACCESS.json")["CFD_access_count"] == 0
    assert d["quality"]["max_accepted_penetration_m"] == d["quality"]["maximum_wall_correction_m"] == 0
    checks["actual_1ms_OFF_zero_drift_no_CFD_original_quality"] = True
    for prefix in ["", "history/pre_top_port/"]:
        base = BASE / prefix
        for rel, record in read(base / "REMOTE_MANIFEST.json")["files"].items():
            p = base / rel
            assert p.stat().st_size == record["bytes"] and sha(p) == record["sha256"], p
    checks["all_remote_files_and_prior_revision_hash_verified"] = True
    prov = read(BASE / "figures/FIGURE_PROVENANCE.json")
    review = read(BASE / "figures/VISUAL_REVIEW.json")
    assert len(prov["figures"]) == len(review["figures"]) == 4 and review["status"] == "PASS"
    assert prov["source_script_sha256"] == sha(ROOT / "scripts/straight_tube_exp1/f10m/figures.py")
    for fig in prov["figures"]:
        for ext, digest in fig["file_sha256"].items():
            assert sha(BASE / "figures" / (fig["figure"] + "." + ext)) == digest
        for rel, digest in fig["source_sha256"].items():
            assert sha(ROOT / rel) == digest
    checks["four_actual_PNG_PDF_pairs_reviewed"] = True
    pg = read(BASE / "progress_document/CHECK.json")
    before = (BASE / "progress_document/before.md").read_bytes()
    after = (BASE / "progress_document/after.md").read_bytes()
    assert sha(pg["backup"]) == pg["before_sha256"] and sha(pg["path"]) == pg["after_sha256"]
    assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
    begin = after.index("\n**F10M ".encode())
    end = after.index("\n---\n\n### 尾声".encode(), begin)
    assert after[:begin] + after[end:] == before
    checks["one_progress_paragraph_backup_old_bytes_retained"] = True
    gate = read(BASE / "analysis/FINAL_GATE.json")
    assert (
        gate["engineering_status"] == "PASS"
        and not gate["training_ready"]
        and not gate["experiment_matched"]
        and gate["three_basis_count"] == 0
    )
    checks["science_limitations_and_NOT_RUN_preserved"] = True
    manifest = BASE / "ARTIFACT_MANIFEST.json"
    if manifest.exists():
        for rel, record in read(manifest)["files"].items():
            assert sha(ROOT / rel) == record["sha256"], rel
        checks["complete_local_artifact_manifest"] = True
    return dict(status="PASS", checks=checks)


def manifest():
    files = [
        p
        for p in BASE.rglob("*")
        if p.is_file() and p.name not in ["ARTIFACT_MANIFEST.json", "VERIFICATION.json"]
    ]
    files += list(CONFIG.glob("*.yaml")) + list((ROOT / "scripts/straight_tube_exp1/f10m").glob("*.py"))
    files += [
        ROOT / p
        for p in [
            "src/mevo_mb/coupling/pdms_block.py",
            "src/mevo_mb/acoustics/spatial_medium.py",
            "src/mevo_mb/acoustics/kwave_backend.py",
            "src/mevo_mb/fast_sim/pdms_engine.py",
            "docs/FAST_PZT_F10M_DEVELOPMENT.md",
            "docs/PDMS_BLOCK_EXPERIMENT_GEOMETRY.md",
            "docs/PDMS_ACOUSTIC_MODEL.md",
        ]
    ]
    data = {str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(set(files))}
    write(BASE / "ARTIFACT_MANIFEST.json", dict(files=data, bytes=sum(v["bytes"] for v in data.values())))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--record", action="store_true")
    a = p.parse_args()
    if a.record:
        manifest()
    result = run()
    if a.record:
        write(BASE / "VERIFICATION.json", result)
    print(result)

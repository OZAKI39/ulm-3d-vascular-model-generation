"""Independently verify similarity glyphs, input identity and complete CFR videos."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path
import subprocess

import h5py
import imageio_ffmpeg
import numpy as np
import pyvista as pv
from scipy.spatial.distance import squareform, pdist
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4l.render_motion import sha, write


def verify(out, decode=True):
    manifest = json.loads((out / "RENDER_MANIFEST.json").read_text())
    assert manifest["display_similarity_scale"] == 24 and not manifest["scientific_data_modified"]
    for p, h in (manifest["protected_inputs"] | manifest["source_sha256"]).items():
        assert sha(p) == h, p
    checks = json.loads((out / "DISPLAY_CHECKS.json").read_text())
    rows = list(csv.DictReader((out / "frame_map.csv").open()))
    report = {}
    for case, video in manifest["videos"].items():
        # NPZ arrays are lazily decompressed on every access. Materialize once
        # before the per-group audit, preserving the exact recorded arrays.
        with np.load(out / f"{case}_RECORDED_SPHERES.npz") as saved:
            native = {key: saved[key] for key in saved.files}
        with np.load(out / f"{case}_DISPLAY_SPHERES.npz") as saved:
            glyph = {key: saved[key] for key in saved.files}
        with h5py.File(video["source"]["trajectory"]) as h:
            for key in ["positions_m", "radii_m", "active", "time_s", "particle_id"]:
                assert np.array_equal(native[key], h[key][:]), key
        for key in ["active", "time_s", "particle_id"]:
            assert np.array_equal(native[key], glyph[key])
        np.testing.assert_array_equal(glyph["radii_m"], native["radii_m"] * 24)
        with gzip.open(out / f"{case}_GROUP_TRANSFORMS.json.gz", "rt") as stream:
            groups = json.load(stream)
        with gzip.open(out / f"{case}_SPHERE_PROJECTIONS.json.gz", "rt") as stream:
            projections = json.load(stream)
        indices = {int(pid): i for i, pid in enumerate(native["particle_id"])}
        max_vector_error = 0.0
        min_radius = np.inf
        counts = []
        intergroup_gap = np.inf
        for k in range(301):
            allids = []
            labels = np.full(len(native["particle_id"]), -1)
            for label, group in enumerate(groups[k]):
                ix = np.array([indices[pid] for pid in group["particle_ids"]])
                labels[ix] = label
                before = native["positions_m"][k, ix]
                after = glyph["positions_m"][k, ix]
                error = np.max(abs((after - after[:1]) - 24 * (before - before[:1])))
                max_vector_error = max(max_vector_error, float(error))
                assert error < 1e-15
                anchor = indices[group["anchor"]]
                expected = (
                    native["positions_m"][k, anchor]
                    + 24 * (before - native["positions_m"][k, anchor])
                    + group["shift_m"]
                )
                np.testing.assert_allclose(after, expected, rtol=0, atol=1e-17)
                assert group["count"] == len(ix) and group["scale"] == 24
                allids.extend(group["particle_ids"])
            active = native["active"][k]
            ids = native["particle_id"][active].tolist()
            assert sorted(allids) == sorted(ids) and len(set(allids)) == len(ids)
            separation = squareform(pdist(glyph["positions_m"][k, active]))
            r = glyph["radii_m"][active]
            gaps = separation - r[:, None] - r[None]
            distinct = labels[active, None] != labels[None, active]
            intergroup_gap = min(intergroup_gap, float(gaps[distinct].min()))
            assert gaps[distinct].min() >= 5e-6 - 1e-12
            counts.append(len(ids))
            if str(k) in projections:
                p = projections[str(k)]
                assert p["submitted_ids"] == ids
                matrix = np.array(p["projection_matrix_px_per_m"])
                offset = np.array(p["projection_offset_px"])
                np.testing.assert_allclose(
                    p["centers_px"], glyph["positions_m"][k, active] @ matrix.T + offset, atol=1e-9, rtol=0
                )
                singular = np.linalg.svd(matrix, compute_uv=False)
                np.testing.assert_allclose(singular, p["pixels_per_m"], rtol=1e-12)
                np.testing.assert_allclose(
                    p["radii_px"], glyph["radii_m"][active] * p["pixels_per_m"], rtol=1e-12
                )
                assert set(p["covered_ids"]).issubset(ids) and len(p["covered_ids"]) > 0
                assert not p["minimum_glyph_size"]
                min_radius = min(min_radius, min(p["radii_px"]))
        identity = Path(video["source"]["trajectory"]).parent.parent / "RUN_IDENTITY.json"
        cfg = json.loads(identity.read_text())["config"]
        geom = BifurcationGeometry(cfg["geometry"])
        surface = next(
            Path(p) for p in manifest["protected_inputs"] if p.endswith("mesh-complete.exterior.vtp")
        )
        geom.bind_mesh_surface(pv.read(surface))
        active = glyph["active"]
        radii = np.broadcast_to(glyph["radii_m"], active.shape)
        gap, _, valid = geom.wall_gap_and_normal(glyph["positions_m"][active], radii[active])
        assert valid.all() and gap.min() >= 0
        scene = checks[case]["scene"]
        assert scene["full_wall_triangles"] == 13978 and scene["open_boundary_edges"] == 90
        assert scene["pzt_body_intersections"] == 0 and scene["reverse_actor_max_rgb_delta"] <= 2
        assert scene["depth_peeling_used"] and scene["back_specular"] == 0.88
        assert len(scene["outward_reflection_lights"]) == 2
        evidence = dict(
            all_301_states_all_active_ids=True,
            unchanged_native_arrays=True,
            maximum_similarity_vector_error_m=max_vector_error,
            minimum_display_wall_gap_m=float(gap.min()),
            minimum_intergroup_display_gap_m=intergroup_gap,
            minimum_display_sphere_diameter_px=2 * min_radius,
            active_counts_first_last=[counts[0], counts[-1]],
            all_spheres_retained=True,
            pzt_geometry_and_transparency_pass=True,
        )
        if manifest["status"] != "PREVIEW":
            assert len(projections) == 301
            timeline = [r for r in rows if r["case"] == case]
            assert len(timeline) == 1201
            for i, row in enumerate(timeline):
                assert int(row["frame"]) == i and int(row["source_index"]) == min(i // 4, 300)
                assert float(row["video_time_s"]) == i / 60
            path = out / video["path"]
            assert sha(path) == video["sha256"]
            reader = imageio_ffmpeg.read_frames(str(path))
            metadata = next(reader)
            reader.close()
            assert metadata["size"] == (3840, 2160) and metadata["fps"] == 60
            decoded = subprocess.run(
                [
                    imageio_ffmpeg.get_ffmpeg_exe(),
                    "-v",
                    "info",
                    "-xerror",
                    "-i",
                    str(path),
                    "-vf",
                    "showinfo",
                    "-an",
                    "-f",
                    "null",
                    "-",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            matches = re.findall(
                r"\[Parsed_showinfo[^\]]*\] n:\s*(\d+) pts:\s*(-?\d+) pts_time:", decoded.stderr
            )
            ticks = re.search(r"config in time_base: (\d+)/(\d+), frame_rate: (\d+)/(\d+)", decoded.stderr)
            assert ticks and ticks.groups()[2:] == ("60", "1")
            timebase = int(ticks[1]) / int(ticks[2])
            assert len(matches) == 1201
            np.testing.assert_array_equal([int(m[0]) for m in matches], np.arange(1201))
            pts = np.array([int(m[1]) for m in matches]) * timebase
            np.testing.assert_allclose(np.diff(pts), 1 / 60, rtol=0, atol=1e-12)
            evidence.update(
                encoding=dict(width=3840, height=2160, fps=60, frames=len(matches), duration_s=1201 / 60),
                constant_frame_spacing=True,
                complete_decode=True,
                sha256=sha(path),
            )
        report[case] = evidence
    return dict(status="PASS", checks=report, raw_scientific_data_modified=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    report = verify(args.output.resolve())
    write(args.output / "DELIVERY_VERIFICATION.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

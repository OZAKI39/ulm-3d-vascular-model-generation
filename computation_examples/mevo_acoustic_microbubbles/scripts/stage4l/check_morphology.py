"""Independent identity, projection, vector geometry and full-video delivery checks."""

from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

import h5py
import imageio.v2 as imageio
import imageio_ffmpeg
import numpy as np
from PIL import Image
import pyvista as pv
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import squareform, pdist

from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4l.render_motion import resolve_input, sha, write


def check(out):
    manifest = json.loads((out / "RENDER_MANIFEST.json").read_text())
    assert not manifest["camera_moves"] and not manifest["particle_interpolation"]
    assert manifest["particle_radius_scale"] == 1 and not manifest["cluster_geometry_transform"]
    for p, digest in manifest["protected_inputs"].items():
        assert sha(p) == digest, p
    rows = list(csv.DictReader((out / "frame_map.csv").open()))
    report = dict(status="PASS", scientific_data_modified=False, videos={})
    for case, video in manifest["videos"].items():
        with h5py.File(video["source"]["trajectory"]) as h:
            original = {k: h[k][:] for k in ["time_s", "positions_m", "radii_m", "active", "particle_id"]}
        with np.load(out / f"{case}_RECORDED_SPHERES.npz") as data:
            for key, values in original.items():
                assert np.array_equal(data[key], values), key
        identity = json.loads((Path(video["source"]["run"]) / "RUN_IDENTITY.json").read_text())
        config = identity["config"]
        surface = config["geometry"]["mesh_surface"]
        mesh = resolve_input(surface, identity["inputs_sha256"][surface])
        geometry = BifurcationGeometry(config["geometry"])
        geometry.bind_mesh_surface(pv.read(mesh))
        x, active, ids, radii = [original[k] for k in ["positions_m", "active", "particle_id", "radii_m"]]
        radii_all = np.broadcast_to(radii, active.shape)
        gap, _, valid = geometry.wall_gap_and_normal(x[active], radii_all[active])
        assert valid.all() and gap.min() >= -1.001e-12
        evidence = json.loads((out / f"{case}_SPHERE_PROJECTIONS.json").read_text())
        assert len(evidence) == len(x) == 301
        first = evidence["0"]
        matrix = np.asarray(first["projection_matrix_px_per_m"])
        offset = np.asarray(first["projection_offset_px"])
        scales = np.linalg.svd(matrix, compute_uv=False)
        np.testing.assert_allclose(scales, scales.mean(), rtol=1e-10)
        max_error = 0.0
        for k in range(301):
            e = evidence[str(k)]
            a = active[k]
            assert e["submitted_ids"] == ids[a].tolist()
            assert set(e["covered_ids"]).issubset(set(ids[a]))
            assert not e["minimum_glyph_size"]
            np.testing.assert_allclose(e["projection_matrix_px_per_m"], matrix, atol=1e-7, rtol=1e-10)
            np.testing.assert_allclose(e["projection_offset_px"], offset, atol=1e-7, rtol=0)
            expected = x[k, a] @ matrix.T + offset
            error = float(np.max(np.abs(np.asarray(e["centers_px"]) - expected)))
            assert error < 1e-7
            max_error = max(max_error, error)
            np.testing.assert_allclose(e["radii_px"], radii[a] * scales.mean(), rtol=1e-10)
            assert (expected > [30, 60]).all() and (expected < [3810, 2100]).all()
            _, labels = connected_components(
                squareform(pdist(x[k, a])) <= config["metrics"]["cluster_threshold_m"], directed=False
            )
            groups = [ids[a][labels == label].tolist() for label in np.unique(labels)]
            groups.sort(key=lambda group: (-len(group), min(group)))
            assert e["largest_cluster_ids"] == groups[0]
        for k in [0, 150, 300]:
            root = ET.parse(out / "vectors" / f"{case}_{k:03d}.svg").getroot()
            circles = root.findall(".//{http://www.w3.org/2000/svg}circle[@data-particle-id]")
            assert sorted(int(c.attrib["data-particle-id"]) for c in circles) == sorted(
                ids[active[k]].tolist()
            )
            mapping = {int(pid): j for j, pid in enumerate(ids)}
            for c in circles:
                i = mapping[int(c.attrib["data-particle-id"])]
                assert np.array_equal(json.loads(c.attrib["data-position-m"]), x[k, i])
                assert float(c.attrib["data-radius-m"]) == radii[i]
                assert abs(float(c.attrib["r"]) - radii[i] * scales.mean()) < 1e-10
                np.testing.assert_allclose(
                    [float(c.attrib["cx"]), float(c.attrib["cy"])],
                    x[k, i] @ matrix.T + offset,
                    atol=1e-7,
                    rtol=0,
                )
        selected = [row for row in rows if row["case"] == case]
        assert len(selected) == 1201
        for f, row in enumerate(selected):
            assert int(row["frame"]) == f and int(row["source_index"]) == min(f // 4, 300)
            assert float(row["source_time_s"]) == original["time_s"][min(f // 4, 300)]
            assert abs(float(row["video_time_s"]) - f / 60) < 1e-13
            assert abs(float(row["requested_time_s"]) - f * 0.005) < 1e-13
        path = out / video["path"]
        assert sha(path) == video["sha256"]
        result = subprocess.run(
            [
                imageio_ffmpeg.get_ffmpeg_exe(),
                "-v",
                "info",
                "-threads",
                "2",
                "-i",
                str(path),
                "-vf",
                "showinfo",
                "-fps_mode",
                "passthrough",
                "-f",
                "null",
                "-",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        assert "error" not in result.stderr.lower(), result.stderr[-2000:]
        timing = re.findall(r"n:\s*(\d+)\s+pts:\s*(-?\d+)\s+pts_time:", result.stderr)
        assert len(timing) == 1201
        frames, pts = np.array(timing, dtype=np.int64).T
        assert np.array_equal(frames, np.arange(1201)) and pts[0] == 0
        tb = re.search(r"config in time_base:\s*(\d+)/(\d+)", result.stderr)
        np.testing.assert_allclose(np.diff(pts) * int(tb[1]) / int(tb[2]), 1 / 60, rtol=0, atol=1e-14)
        reader = imageio.get_reader(path)
        meta = reader.get_meta_data()
        assert meta["fps"] == 60 and list(meta["size"]) == [3840, 2160]
        errors = {}
        for frame in [0, 600, 1200]:
            decoded = reader.get_data(frame)
            reference = np.asarray(Image.open(out / "keyframes" / f"{case}_{frame:06d}.png"))
            mae = float(np.abs(decoded.astype(float) - reference).mean())
            assert mae < 1.5
            errors[str(frame)] = mae
        reader.close()
        report["videos"][case] = dict(
            frames=1201,
            fps=60,
            resolution=[3840, 2160],
            full_decode="PASS",
            timestamps="Exact 1/60 s; no missing encoded frames",
            stored_states=301,
            every_active_id_submitted=True,
            saved_positions_and_radii_bitwise_exact=True,
            camera_fixed=True,
            common_isotropic_physical_scale=True,
            projected_position_max_error_px=max_error,
            all_centers_inside=True,
            minimum_wall_gap_m=float(gap.min()),
            vector_identity_and_geometry="PASS",
            glyph_diameter_px=float(2 * radii[0] * scales.mean()),
            snapshot_MAE=errors,
        )
        print(json.dumps(dict(case=case, status="PASS", frames=1201)), flush=True)
    write(out / "DELIVERY_VERIFICATION.json", report)
    manifest["status"] = "VERIFIED"
    write(out / "RENDER_MANIFEST.json", manifest)
    write(
        out / "DELIVERY_SHA256.json",
        {
            str(p.relative_to(out)): sha(p)
            for p in sorted(out.rglob("*"))
            if p.is_file() and p.name not in ["DELIVERY_SHA256.json", "render.log", "check.log"]
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    check(parser.parse_args().output.resolve())

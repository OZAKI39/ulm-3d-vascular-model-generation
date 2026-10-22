"""Verify video timing, decoded images, recorded populations and source integrity."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

import h5py
import imageio.v2 as imageio
import imageio_ffmpeg
import numpy as np
from PIL import Image


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def check_physical_centers(out, key, video, display_times):
    import pyvista as pv
    from mevo_mb.coupling.bifurcation import BifurcationGeometry
    from scripts.stage4l.render_motion import resolve_input, style

    run = Path(video["source"]["run"])
    identity = json.loads((run / "RUN_IDENTITY.json").read_text())
    cfg = identity["config"]
    with h5py.File(video["source"]["trajectory"]) as h:
        times, x, ids, radii, active = [
            h[k][:] for k in ["time_s", "positions_m", "particle_id", "radii_m", "active"]
        ]
    with np.load(out / video["display_samples"]) as s:
        q, drawn, visible, correction = [
            s[k] for k in ["time_s", "positions_m", "active", "wall_adjustment_m"]
        ]
        assert np.array_equal(s["particle_id"], ids)
    assert np.array_equal(q, display_times)
    assert np.array_equal(drawn[::4], x), "Saved particle centers were displaced"
    assert np.array_equal(visible[::4], active)
    assert not ((~visible[:-1]) & visible[1:]).any()
    branches = json.loads((run / "BRANCH_METRICS.json").read_text())
    mapping = {int(pid): i for i, pid in enumerate(ids)}
    exits = np.full(len(ids), np.inf)
    for e in branches["crossings"]:
        exits[mapping[e["particle_id"]]] = e["crossing_time_s"]
    assert np.array_equal(visible, q[:, None] < exits)
    straight = style.base.sample_case(dict(t=times, x=x, id_index=mapping, events=branches["crossings"]), q)[
        0
    ]
    straight[::4] = x
    actual_correction = np.linalg.norm(drawn - straight, axis=2)
    np.testing.assert_allclose(actual_correction, correction, atol=1e-18, rtol=0)
    assert not correction[::4].any() and correction.max() < 10e-6
    surface = cfg["geometry"]["mesh_surface"]
    mesh = resolve_input(surface, identity["inputs_sha256"][surface])
    geometry = BifurcationGeometry(cfg["geometry"])
    geometry.bind_mesh_surface(pv.read(mesh))
    r = np.broadcast_to(radii, visible.shape)[visible]
    gap, _, valid = geometry.wall_gap_and_normal(drawn[visible], r)
    assert valid.all() and np.all(gap + r > 0) and gap.min() >= -1.001e-12
    checks = json.loads((out / "DISPLAY_GEOMETRY_CHECKS.json").read_text())[key]
    assert checks["frames_checked_for_glyph_clipping"] == video["frames"]
    assert checks["outside_lumen_glyph_pixels"] == 0
    common = visible[1:] & visible[:-1]
    maximum_step = float(np.linalg.norm(np.diff(drawn, axis=0), axis=2)[common].max())
    return dict(
        saved_positions_bitwise_exact=True,
        independently_checked_active_samples=int(visible.sum()),
        active_centers_outside_count=int((gap + r < 0).sum()),
        minimum_surface_gap_m=float(gap.min()),
        max_display_wall_adjustment_m=float(correction.max()),
        max_display_frame_displacement_m=maximum_step,
        group_coordinate_magnification=False,
        raster_cross_dissolve=False,
        glyph_clipping_all_frames=checks["frames_checked_for_glyph_clipping"],
    )


def check(out):
    manifest = json.loads((out / "RENDER_MANIFEST.json").read_text())
    physical_centers = manifest.get("display_geometry") == "physical_centers"
    for p, digest in (manifest["protected_inputs"] | manifest["renderer_sources"]).items():
        assert sha(p) == digest, p
    rows = list(csv.DictReader((out / "frame_map.csv").open()))
    states = json.loads((out / "STATE_GROUPS.json").read_text())
    report = dict(
        status="PASS",
        videos={},
        source_hashes_unchanged=True,
        scope="Uniform encoding and display clock do not certify physical accuracy. Physically stationary/subpixel motion is not a dropped frame.",
    )
    for key, video in manifest["videos"].items():
        path = out / video["path"]
        assert sha(path) == video["sha256"]
        n = video["frames"]
        selected = [row for row in rows if row["video"] == key]
        assert len(selected) == n
        display_times = np.array([float(row["display_time_s"]) for row in selected])
        video_times = np.array([float(row["video_time_s"]) for row in selected])
        expected_step = (video["source"]["time_end_s"] - video["source"]["time_start_s"]) / (n - 1)
        np.testing.assert_allclose(np.diff(display_times), expected_step, rtol=0, atol=2e-15)
        np.testing.assert_allclose(video_times, np.arange(n) / 60, atol=1e-14, rtol=0)
        assert display_times[0] == video["source"]["time_start_s"]
        assert abs(display_times[-1] - video["source"]["time_end_s"]) < 1e-12
        with h5py.File(video["source"]["trajectory"], "r") as f:
            times, active, ids = f["time_s"][:], f["active"][:], f["particle_id"][:]
        for frame, row in enumerate(selected):
            assert int(row["frame"]) == frame
            left, right = int(row["source_left"]), int(row["source_right"])
            alpha = float(row["image_blend_fraction"])
            assert float(row["source_left_time_s"]) == times[left]
            assert float(row["source_right_time_s"]) == times[right]
            if video["kind"] == "dense_probe":
                assert left == right and alpha == 0
                assert left == np.searchsorted(times, display_times[frame], side="right") - 1
            else:
                assert left == frame // 4 and right == min(left + 1, 300)
                if physical_centers:
                    assert alpha == 0
                    coordinate_fraction = float(row["coordinate_fraction"])
                    assert coordinate_fraction == (frame % 4) / 4
                else:
                    assert alpha == (frame % 4) / 4
                    coordinate_fraction = alpha
                assert (
                    abs(
                        display_times[frame]
                        - ((1 - coordinate_fraction) * times[left] + coordinate_fraction * times[right])
                    )
                    < 1e-12
                )
        for state in states[key]:
            source = state["source_state"]
            actual_ids = sorted(ids[active[source]].tolist())
            assert sorted(state["visible_particle_ids"]) == actual_ids
            group_ids = [pid for group in state["groups"] for pid in group["particle_ids"]]
            assert len(group_ids) == len(set(group_ids))
            assert sorted(group_ids) == actual_ids
            assert max(g["count"] for g in state["groups"]) == state["stats"]["largest_cluster_size"]
            assert all(g["count"] == len(g["particle_ids"]) for g in state["groups"])
        # Decode without frame-rate repair, and inspect every real presentation timestamp.
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
        assert len(timing) == n
        frame_numbers, pts = np.array(timing, dtype=np.int64).T
        assert np.array_equal(frame_numbers, np.arange(n)) and pts[0] == 0
        tb = re.search(r"config in time_base:\s*(\d+)/(\d+)", result.stderr)
        tick = int(tb[1]) / int(tb[2])
        np.testing.assert_allclose(np.diff(pts) * tick, 1 / 60, atol=1e-14, rtol=0)
        expected = {int(p.stem.rsplit("_", 1)[1]): p for p in (out / "keyframes").glob(f"{key}_*.png")}
        reader = imageio.get_reader(path)
        meta = reader.get_meta_data()
        assert meta["fps"] == 60 and list(meta["size"]) == [1920, 1080]
        snapshot_errors = {}
        for frame, p in sorted(expected.items()):
            im = reader.get_data(frame)
            reference = np.asarray(Image.open(p))
            mae = float(np.abs(im.astype(float) - reference).mean())
            assert mae < 2, (key, frame, mae)
            assert im[1010:].max() < 12, "Unwanted footer"
            assert im[220:267, 60:450].max() > 150, "Missing upper-left clock"
            assert im[15:80, 1600:1910].max() < 12, "Unwanted top-right clock"
            snapshot_errors[str(frame)] = mae
        reader.close()
        report["videos"][key] = dict(
            frames=n,
            fps=60,
            resolution=[1920, 1080],
            timestamps="uniform 1/60 s, no missing encoded frames",
            full_decode="PASS",
            snapshot_MAE_0_to_255=snapshot_errors,
            constant_physical_time_step_s=expected_step,
            playback_speed=video["playback_speed"],
            exact_active_ID_partition=True,
            source_mapping="PASS",
            kind=video["kind"],
        )
        if physical_centers:
            report["videos"][key]["display_geometry"] = check_physical_centers(out, key, video, display_times)
        print(json.dumps(dict(video=key, status="PASS", frames=n)), flush=True)
    assert len(rows) == sum(v["frames"] for v in manifest["videos"].values())
    (out / "DELIVERY_VERIFICATION.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest["status"] = "VERIFIED"
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    artifacts = {
        str(p.relative_to(out)): sha(p)
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name not in ["DELIVERY_SHA256.json", "render.log", "check.log"]
    }
    (out / "DELIVERY_SHA256.json").write_text(json.dumps(artifacts, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    check(parser.parse_args().output.resolve())

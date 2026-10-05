"""Decode and verify Stage 4 overview/zoom deliveries against their saved-state images."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

import imageio.v2 as imageio
import imageio_ffmpeg
import numpy as np
from PIL import Image


def sha(path):
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def check(out):
    manifest = json.loads((out / "RENDER_MANIFEST.json").read_text())
    single_view = manifest.get("layout") == "single_overview"
    for path, digest in (manifest["protected_inputs"] | manifest["renderer_sources"]).items():
        assert sha(path) == digest, path
    rows = list(csv.DictReader((out / "frame_map.csv").open()))
    exact = manifest["temporal_display"]["exact_state_video_frames"]
    assert len(exact) == 301 and exact[0] == 0 and exact[-1] == 1400
    assert np.array_equal(np.diff(exact), [24] * 10 + [4] * 290)
    report = dict(status="PASS", video_checks={}, frame_map_rows=len(rows))
    for name, info in manifest["videos"].items():
        path = out / info["path"]
        assert sha(path) == info["sha256"]
        # showinfo exposes actual presentation timestamps; passthrough avoids FPS repair.
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
        assert len(timing) == 1401, len(timing)
        frame_numbers, pts = np.array(timing, dtype=np.int64).T
        assert np.array_equal(frame_numbers, np.arange(1401))
        tb = re.search(r"config in time_base:\s*(\d+)/(\d+)", result.stderr)
        tick_s = int(tb[1]) / int(tb[2])
        assert pts[0] == 0 and np.allclose(np.diff(pts) * tick_s, 1 / 60, rtol=0, atol=1e-14)
        reader = imageio.get_reader(path)
        meta = reader.get_meta_data()
        assert meta["fps"] == 60 and list(meta["size"]) == info["resolution"]
        # Decode all frames again for content checks: clock excluded from motion checksum.
        last_hash, consecutive, longest, repeated, last_image = None, 0, 0, 0, None
        snapshot_errors, decoded = {}, 0
        expected = {}
        if name in ["B2", "B3"]:
            for png in (out / "keyframes").glob(f"{name}_*.png"):
                suffix = png.stem.removeprefix(name + "_")
                frame = int(suffix.split("_")[-1]) if suffix.startswith("transition_") else exact[int(suffix)]
                expected[frame] = png
        for frame, im in enumerate(reader):
            decoded += 1
            scene = im[240:930] if name in ["B2", "B3"] else im[120:465]
            digest = hashlib.blake2b(scene.tobytes(), digest_size=16).digest()
            if digest == last_hash:
                repeated += 1
                consecutive += 1
            else:
                consecutive = 0
            longest = max(longest, consecutive)
            last_hash, last_image = digest, im
            if frame in expected:
                reference = np.asarray(Image.open(expected[frame]))
                mae = float(np.abs(im.astype(float) - reference).mean())
                assert mae < 2, (name, frame, mae)
                snapshot_errors[str(frame)] = mae
            if name == "B3" and not single_view:
                r = rows[1401 + frame]
                if r["anchor_active"] == "False":
                    # Empty zoom must not retain blue/mint bubble glyphs after crossing.
                    crop = im[400:750, 1400:1660].astype(float)
                    blue = (crop[:, :, 2] > 90) & (crop[:, :, 2] > crop[:, :, 0] * 1.8) & (crop[:, :, 1] > 45)
                    assert blue.sum() == 0, (frame, int(blue.sum()))
        reader.close()
        assert decoded == 1401
        assert longest < 3, (name, "Repeated scene frames", longest)
        if name in ["B2", "B3"]:
            assert last_image[1010:].max() < 12, "Unwanted footer content"
        report["video_checks"][name] = dict(
            frames=decoded,
            fps=60,
            timestamps="uniform, no gaps",
            full_decode="PASS",
            resolution=info["resolution"],
            repeated_scene_frames=repeated,
            longest_repeat_run=longest,
            snapshot_MAE_0_to_255=snapshot_errors,
        )
        print(json.dumps(dict(case=name, decoded_frames=decoded, status="PASS")), flush=True)
    assert len(rows) == 2802
    for name in ["B2", "B3"]:
        selected = [r for r in rows if r["case"] == name]
        assert len(selected) == 1401
        times = np.array([float(r["display_time_s"]) for r in selected])
        assert times[0] == 0 and times[-1] == 6 and np.all(np.diff(times) > 0)
        for frame, r in enumerate(selected):
            left, right = int(r["source_left"]), int(r["source_right"])
            alpha = float(r["image_blend_fraction"])
            assert int(r["frame"]) == frame and float(r["video_time_s"]) == frame / 60
            assert 0 <= alpha < 1 and right in [left, left + 1]
            assert abs(float(r["display_time_s"]) - 0.02 * (left + alpha)) < 1e-12
        for state, frame in enumerate(exact):
            assert int(selected[frame]["source_left"]) == state
            assert float(selected[frame]["image_blend_fraction"]) == 0
        stats = manifest["states"][name]
        assert len(stats) == 301
        assert all(s["A"] + s["B"] + s["active_particles"] == 256 for s in stats)
        assert all(set(s["cluster_particle_ids"]) <= set(s["visible_particle_ids"]) for s in stats)
        final = stats[-1]
        assert [final["A"], final["B"], final["active_particles"]] == (
            [76, 0, 180] if name == "B2" else [0, 88, 168]
        )
        if single_view:
            import h5py

            original_file = next(
                Path(path)
                for path in manifest["protected_inputs"]
                if path.endswith("particles/trajectories.h5") and f"_{name}_seed" in path
            )
            with h5py.File(original_file) as source:
                positions, ids, active = (source[key][:] for key in ["positions_m", "particle_id", "active"])
            mapping = {int(pid): i for i, pid in enumerate(ids)}
            for k, state in enumerate(stats):
                groups = state["overview_groups"]
                members = [pid for group in groups for pid in group["particle_ids"]]
                assert len(members) == len(set(members)) == state["active_particles"]
                assert set(members) == set(ids[active[k]])
                assert len(groups) == state["stats"]["cluster_count"]
                assert max(g["count"] for g in groups) == state["stats"]["largest_cluster_size"]
                for group in groups:
                    indices = [mapping[pid] for pid in group["particle_ids"]]
                    assert group["count"] == len(indices)
                    np.testing.assert_allclose(
                        group["center_m"], positions[k, indices].mean(axis=0), rtol=0, atol=1e-15
                    )
                    assert group["marker_area_pt2"] == 18 * len(indices)
            layout = manifest["screen_checks"][name]
            assert layout["three_dimensional_axes"] == 1
            assert not layout["zoom_panel"] and not layout["roi_connectors"]
            assert layout["scale_relative_to_previous_left_view"] > 1.2
            report.setdefault("single_view_checks", {})[name] = dict(
                frames_checked=301,
                exact_active_partition=True,
                true_group_centers=True,
                marker_area_proportional_to_count=True,
                one_scene_only=True,
                enlargement=layout["scale_relative_to_previous_left_view"],
            )
    html = (out / "OPEN_RESULTS.html").read_text()
    for link in re.findall(r'(?:href|src|poster)="([^"]+)"', html):
        if link == "DELIVERY_VERIFICATION.json":
            continue  # This check creates the linked report below, after validation.
        assert (out / link).is_file(), link
    report["scientific_inputs"] = (
        "All original hashes unchanged; all 301 states retain active population and recorded groups"
        if single_view
        else "All original hashes unchanged; all 301 states retain full tracked component"
    )
    report["source_state_semantics"] = (
        "Cross-dissolve images only; actual saved coordinates and radii unchanged"
    )
    (out / "DELIVERY_VERIFICATION.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest["status"] = "PASS"
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.copy2(__file__, out / "provenance" / Path(__file__).name)
    hashes = {
        str(p.relative_to(out)): dict(bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name not in ["DELIVERY_SHA256.json", "render.log"]
    }
    (out / "DELIVERY_SHA256.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    check(parser.parse_args().output.resolve())

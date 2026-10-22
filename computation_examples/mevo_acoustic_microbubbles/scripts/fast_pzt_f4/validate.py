"""Validate closed science arrays, compact copies and full MP4 decoding/time stamps."""

import hashlib
import json
import re
import subprocess
import numpy as np
import imageio_ffmpeg
from .common import BASE, ROOT, write, js
from .analyze import B_NAMES, C_NAMES, read, R


def main():
    bounds = np.array(js(BASE / "GEOMETRY_FIELD_OVERLAP.json")["usable_center_bounds_m"])
    rows = {}
    for group, names in [("aggregation_v1", B_NAMES), ("controls_v1", C_NAMES)]:
        for name in names:
            path = BASE / group / name
            d = read(path)
            x = d["positions_m"][:, 0]
            assert np.isfinite(x).all() and d["active"].all()
            assert (x >= bounds[:, 0]).all() and (x <= bounds[:, 1]).all()
            np.testing.assert_array_equal(d["particle_id"], np.arange(32))
            assert np.all(d["radii_m"] == R)
            i, j = np.triu_indices(32, 1)
            gap = np.linalg.norm(x[:, i] - x[:, j], axis=-1) - 2 * R
            assert gap.min() >= -0.01 * R
            if "wall_on" in name:
                assert d["wall_gap_m"].min() >= -0.01 * R
            assert json.loads(js(path / "SAVED_METADATA.json")["text_datasets"]["boundary_events_json"]) == []
            keys = []
            for copy in path.glob("SAVED_*.npz"):
                with np.load(copy) as data:
                    for key in data.files:
                        np.testing.assert_array_equal(data[key], d[key])
                        keys.append(key)
            rows[group + "/" + name] = dict(
                saved_frames=len(x),
                minimum_saved_pair_gap_m=float(gap.min()),
                compact_array_count=len(keys),
                compact_arrays_bitwise_equal=True,
                finite=True,
                all_centres_in_stored_field=True,
                all_IDs_active=True,
            )
    # All inherited executable files remain unchanged, not merely the acoustic formula names.
    changes = subprocess.check_output(
        [
            "git",
            "diff",
            "--name-status",
            "3eeb94eb885f5df4f0e214ba9b7904dbdecec744",
            "--",
            "src/mevo_mb/fast_sim",
        ],
        cwd=ROOT,
        text=True,
    )
    assert all(line.startswith("A\t") for line in changes.splitlines()), changes
    write(
        BASE / "ARRAY_VALIDATION.json",
        dict(status="PASS", rows=rows, original_fast_sim_files_unchanged=True, git_diff=changes),
    )
    file = BASE / "visualization/F4_Y_vessel_aggregation_and_PZT_control_60fps.mp4"
    cmd = [
        imageio_ffmpeg.get_ffmpeg_exe(),
        "-hide_banner",
        "-i",
        str(file),
        "-vf",
        "showinfo",
        "-f",
        "null",
        "-",
    ]
    run = subprocess.run(cmd, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-3000:]
    frames = re.findall(r" n:\s*(\d+)\s+pts:\s*(-?\d+)\s+pts_time:([-+0-9.e]+)", run.stderr)
    base = re.search(r"config in time_base:\s*(\d+)/(\d+)", run.stderr)
    assert base and len(frames) == 721, len(frames)
    fps = js(BASE / "visualization/VIDEO_IDENTITY.json")["fps"]
    pts = np.array([int(f[1]) for f in frames])
    t = pts * int(base[1]) / int(base[2])
    np.testing.assert_array_equal([int(f[0]) for f in frames], np.arange(721))
    np.testing.assert_allclose(t, np.arange(721) / fps, rtol=0, atol=1e-12)
    assert len(np.unique(np.diff(pts))) == 1
    with np.load(BASE / "visualization/FRAME_MAP.npz") as m:
        np.testing.assert_allclose(m["physical_time_s"], np.arange(721) * 0.001, rtol=0, atol=1e-14)
        assert len(m["video_frame"]) == 721
    write(
        BASE / "visualization/VIDEO_VALIDATION.json",
        dict(
            status="PASS",
            decoded_frames=len(frames),
            fps=fps,
            timestamp_increment_ticks=int(np.diff(pts)[0]),
            time_base=[int(base[1]), int(base[2])],
            duration_s=721 / fps,
            physical_time_range_s=[0, 0.720],
            constant_PTS=True,
            sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
            bytes=file.stat().st_size,
            encoding="H.264/yuv420p; 1920x1080; faststart",
            decode_error=False,
            evidence="Every frame decoded with ffmpeg showinfo; exact integer PTS cadence checked; actual source frame map verified.",
        ),
    )
    print("Arrays and full video decode PASS")


if __name__ == "__main__":
    main()

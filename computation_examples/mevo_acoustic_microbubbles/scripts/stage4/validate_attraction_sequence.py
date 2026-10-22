"""Actual-frame and continuity checks for the sequential attraction delivery."""

from fractions import Fraction
import re
import subprocess

import h5py
import imageio_ffmpeg
import numpy as np
from PIL import Image
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform

from .visual_attraction_sequence import OUT, ROOT, RUN, STEM, inputs, style


def main():
    data, cfg = inputs()
    summary = style.js(RUN / "SUMMARY.json")
    record = summary["switch_records"][0]
    assert record["time_s"] == 0.3 and record["positions_unchanged"]
    assert record["discarded_graphs"] > 0
    assert record["primary_nonzero"] and record["secondary_addition_exactly_zero"]
    assert summary["quality"]["error_bits"] == 0
    for key in ["amplitude_Pa", "phase_rad", "enabled"]:
        assert cfg["events"][0][key] == cfg["events"][1][key]
    with h5py.File(style.F4 / "aggregation_v1/parent_wall_on/trajectory.h5") as h:
        np.testing.assert_array_equal(data["x"][:601:2], h["positions_m"][:1501:5, 0])
    np.testing.assert_allclose(np.diff(data["time"]), 0.0005, rtol=1e-12, atol=1e-16)
    rg = np.sqrt(np.mean(np.sum((data["x"] - data["x"].mean(1)[:, None]) ** 2, axis=2), axis=1))
    np.testing.assert_allclose(rg, data["rg"], rtol=1e-13, atol=1e-19)
    for i, x in enumerate(data["x"]):
        count, labels = connected_components(
            squareform(pdist(x)) - 2 * style.R <= 0.1 * style.R, directed=False
        )
        assert count == data["nc"][i] and np.bincount(labels).max() == data["nmax"][i]
    identity = style.js(RUN / "IDENTITY.json")
    assert identity["runner_sha256"] == style.sha(ROOT / "scripts/stage4/visual_sequential_ablation.py")
    for name, expected in identity["engine_sha256"].items():
        assert style.sha(style.F4.resolve().parents[2] / "src/mevo_mb/fast_sim" / name) == expected

    preserved = 0
    for version in ["f3", "f4"]:
        p = ROOT / "results/fast_pzt" / version
        for row in style.js(p / "ARTIFACT_MANIFEST.json")["artifacts"]:
            assert style.sha(p / row["relative_path"]) == row["sha256"]
            preserved += 1
    old = style.OUT
    for row in style.js(old / "ARTIFACT_MANIFEST.json")["files"]:
        assert style.sha(old / row["path"]) == row["sha256"]
        preserved += 1
    result = subprocess.run(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-nostdin",
            "-threads",
            "2",
            "-i",
            str(OUT / (STEM + ".mp4")),
            "-vf",
            "showinfo",
            "-an",
            "-f",
            "null",
            "-",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    (OUT / "decode.log").write_text(result.stderr)
    rows = re.findall(r"n:\s*(\d+)\s+pts:\s*(-?\d+)\s+pts_time:([^\s]+)", result.stderr)
    assert len(rows) == 1201
    np.testing.assert_array_equal([int(r[0]) for r in rows], np.arange(1201))
    pts = np.array([int(r[1]) for r in rows])
    num, den = map(int, re.search(r"config in time_base:\s*(\d+)/(\d+)", result.stderr).groups())
    assert np.unique(np.diff(pts)).size == 1 and pts[0] == 0
    assert Fraction(int(np.diff(pts)[0]) * num, den) == Fraction(1, 60)
    mapping = np.load(OUT / "FRAME_MAP.npz")
    np.testing.assert_array_equal(mapping["time_s"], data["time"])
    np.testing.assert_array_equal(mapping["source_frame"], np.arange(1201))
    np.testing.assert_array_equal(mapping["secondary_enabled"], np.arange(1201) < 600)
    for index in [0, 200, 599, 600, 1200]:
        with Image.open(OUT / (STEM + f"_frame_{index:04d}.png")) as im:
            assert im.size == (1920, 1080)
            assert np.array_equal(np.asarray(im)[0, 0, :3], [0, 0, 0])
    render = style.js(OUT / "RENDER_SOURCE.json")
    assert render["renderer"] == style.sha(ROOT / "scripts/stage4/visual_attraction_sequence.py")
    assert render["style_dependency"] == style.sha(ROOT / "scripts/stage4/visual_mechanism_stories.py")
    report = dict(
        status="PASS",
        data_frames=1201,
        decoded_frames=1201,
        fps="60/1",
        resolution=[1920, 1080],
        time_interval_s=0.0005,
        same_particles=True,
        same_PZT=True,
        positions_preserved_at_switch=True,
        ON_graphs_discarded=True,
        on_prefix_bitwise_reference_matches=301,
        old_files_unchanged=preserved,
        switch_time_s=0.3,
        switch_video_time_s=10.0,
        switch_largest=int(data["nmax"][600]),
        final_largest=int(data["nmax"][-1]),
        switch_rg_um=float(rg[600] * 1e6),
        final_rg_um=float(rg[-1] * 1e6),
        maximum_after_switch_displacement_um=float(
            np.linalg.norm(data["x"][600:] - data["x"][600], axis=2).max() * 1e6
        ),
        video_sha256=style.sha(OUT / (STEM + ".mp4")),
        scope="Software/data/render checks; not independent physical validation",
    )
    style.write(OUT / "VALIDATION.json", report)
    print(report)


if __name__ == "__main__":
    main()

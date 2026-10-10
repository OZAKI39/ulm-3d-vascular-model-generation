"""Validate the actual data, preservation, and complete video frame timing."""

import json
import re
import subprocess
from fractions import Fraction

import h5py
import imageio_ffmpeg
import numpy as np
from PIL import Image
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform

from .visual_mechanism_stories import ROOT, OUT, F4, R, load, sha, js, write


def main():
    paths, data = load()
    report = {"sources": {}, "archive_preservation": {}, "videos": {}, "figures": {}}
    for version in ["f3", "f4"]:
        directory = ROOT / "results/fast_pzt" / version
        artifacts = js(directory / "ARTIFACT_MANIFEST.json")["artifacts"]
        for item in artifacts:
            assert sha(directory / item["relative_path"]) == item["sha256"], item["relative_path"]
        report["archive_preservation"][version] = {"checked": len(artifacts), "mismatches": 0}

    ca = js(paths["aggregation_on"] / "CONFIG.json")
    cb = js(paths["aggregation_off"] / "CONFIG.json")
    for k, v in ca.items():
        if k != "scenario":
            assert cb[k] == v, k
    assert cb["secondary_force_enabled"] is False
    ident = js(paths["aggregation_off"] / "IDENTITY.json")
    assert ident["runner_sha256"] == sha(ROOT / "scripts/stage4/visual_evidence_ablation.py")
    assert ident["matched_secondary_on_sha256"] == sha(paths["aggregation_on"] / "trajectory.h5")
    archived_engine = F4.resolve().parents[2] / "src/mevo_mb/fast_sim"
    for name, digest in ident["engine_sha256"].items():
        assert sha(archived_engine / name) == digest, name
    report["ablation"] = dict(
        same_initial_positions=True,
        same_PZT_drive=True,
        same_dt_and_output_times=True,
        only_secondary_force_addition_disabled=True,
        engine_matches_archived_F4=True,
        runner_matches_executed_source=True,
    )

    for key, d in data.items():
        recomputed_rg = np.sqrt(np.mean(np.sum((d["x"] - d["x"].mean(1)[:, None]) ** 2, axis=2), axis=1))
        np.testing.assert_allclose(recomputed_rg, d["rg"], rtol=1e-13, atol=1e-19)
        for i, positions in enumerate(d["x"]):
            adjacency = squareform(pdist(positions)) - 2 * R <= 0.1 * R
            count, labels = connected_components(adjacency, directed=False)
            assert count == d["nc"][i]
            assert np.bincount(labels).max() == d["nmax"][i]
        summary = js(paths[key] / "SUMMARY.json")
        assert summary["quality"]["error_bits"] == 0
        with h5py.File(paths[key] / "trajectory.h5") as h:
            if "active" in h:
                assert h["active"][:].all()
        report["sources"][key] = dict(
            verified_frames=len(d["time"]),
            all_32_ids_retained=True,
            radius_m=R,
            metrics_recomputed=True,
            quality=summary["numerical_quality"],
            source_sha256=sha(paths[key] / "trajectory.h5"),
        )

    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    for kind, stem, source in [
        ("aggregation", "Video_1_attraction_on_off_60fps", "aggregation_on"),
        ("translation", "Video_2_PZT_translation_shape_limit_60fps", "translation_baseline"),
    ]:
        result = subprocess.run(
            [
                ffmpeg,
                "-nostdin",
                "-threads",
                "2",
                "-i",
                str(OUT / (stem + ".mp4")),
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
        (OUT / (kind + "_decode.log")).write_text(result.stderr)
        rows = re.findall(r"n:\s*(\d+)\s+pts:\s*(-?\d+)\s+pts_time:([^\s]+)", result.stderr)
        assert len(rows) == 601, len(rows)
        np.testing.assert_array_equal([int(r[0]) for r in rows], np.arange(601))
        pts = np.array([int(r[1]) for r in rows])
        num, den = map(int, re.search(r"config in time_base:\s*(\d+)/(\d+)", result.stderr).groups())
        assert np.unique(np.diff(pts)).size == 1 and pts[0] == 0
        assert Fraction(int(np.diff(pts)[0]) * num, den) == Fraction(1, 60)
        assert "s:1920x1080" in result.stderr
        mapping = np.load(OUT / (stem + "_FRAME_MAP.npz"))
        stride = 5 if kind == "aggregation" else 1
        np.testing.assert_array_equal(mapping["source_frame_index"], np.arange(601) * stride)
        np.testing.assert_array_equal(
            mapping["physical_time_s"], data[source]["time"][mapping["source_frame_index"]]
        )
        np.testing.assert_allclose(
            np.diff(mapping["physical_time_s"]), 0.001 if stride == 5 else 0.0002, rtol=1e-12, atol=1e-16
        )
        render = js(OUT / (stem + "_IDENTITY.json"))
        assert not render["position_interpolation"] and not render["spatial_rescaling"]
        assert render["bubble_radius_scale"] == 1
        assert js(OUT / ("RENDER_SOURCE_" + kind + ".json"))["sha256"] == sha(
            ROOT / "scripts/stage4/visual_mechanism_stories.py"
        )
        report["videos"][kind] = dict(
            complete_decode="PASS",
            frames=601,
            size=[1920, 1080],
            fps="60/1",
            uniform_encoded_frame_intervals=True,
            uniform_physical_time_intervals=True,
            duration_s=601 / 60,
            source_frame_map_checked=True,
            interpolation=False,
            spatial_rescaling=False,
            sha256=sha(OUT / (stem + ".mp4")),
        )

    for stem in ["Fig_1_attraction_on_off", "Fig_2_PZT_translation_shape_limit"]:
        with Image.open(OUT / (stem + ".png")) as im:
            assert im.size == (4800, 2550)
            assert np.asarray(im.convert("RGB")).std() > 10
            report["figures"][stem] = dict(
                size=im.size, png_dpi=im.info.get("dpi"), formats=["png", "pdf", "svg"]
            )
        for extension in ["pdf", "svg"]:
            assert (OUT / (stem + "." + extension)).stat().st_size > 1000
    assert js(OUT / "RENDER_SOURCE_figures.json")["sha256"] == sha(
        ROOT / "scripts/stage4/visual_mechanism_stories.py"
    )
    for name in ["visual_evidence_ablation.py", "visual_mechanism_stories.py"]:
        assert sha(ROOT / "scripts/stage4" / name) == sha(OUT / "source" / name)
    report["status"] = "PASS"
    report["scientific_scope"] = (
        "Internal software/data/render validation, not independent validation of near-contact physics."
    )
    write(OUT / "VALIDATION.json", report)
    print(
        json.dumps(
            {
                "status": "PASS",
                "source_frames": sum(v["verified_frames"] for v in report["sources"].values()),
                "decoded_video_frames": 1202,
                "unchanged_archive_files": sum(v["checked"] for v in report["archive_preservation"].values()),
            }
        )
    )


if __name__ == "__main__":
    main()

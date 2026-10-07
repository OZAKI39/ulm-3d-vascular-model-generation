"""Equal-camera old/new +5 ms video, real radii and actual stored frames only."""

import argparse
from pathlib import Path
import imageio.v2 as imageio
import numpy as np
import pyvista as pv
from scripts.stage5a_r3.common import BASE, R2, read, write, sha
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r3.figures import history, BLUE
from scripts.stage5a_r3.video_audit import audit
from scripts.stage5a_r2.motion_visuals import wall_path, glyph, decorate


def main(case):
    ledger = read(BASE / "final/CAMPAIGN_R3_V1/RESULT.json")["ledger"]
    chosen = [
        next(
            r
            for r in ledger
            if r["case"] == case
            and r["duration_s"] == 0.005
            and r["kind"] == kind
            and r["repeat"] == 0
            and r.get("run")
        )
        for kind in ["reference", "nominal"]
    ]
    assert all(r["summary"]["status"] == "COMPLETE" for r in chosen)
    paths = [local(r["run"]) / "trajectories.h5" for r in chosen]
    states = [history(p) for p in paths]
    for key in ["particle_id", "radii_m"]:
        assert np.array_equal(states[0][key], states[1][key])
    assert np.array_equal(states[0]["positions_m"][0], states[1]["positions_m"][0], equal_nan=True)
    assert all(np.all(s["active"] == s["active"][0]) for s in states)
    group = np.load(R2 / "comparisons" / f"{case}_0.001_v4a_legacy_reference.npz")["particle_ids"]
    masks = [s["active"][0] for s in states]
    assert np.array_equal(*masks)
    r = states[0]["radii_m"][masks[0]]
    assert np.all(r == r[0]), "Use a per-ID sphere mapper for unequal radii"
    radius = float(r[0])
    selected = np.isin(states[0]["particle_id"], group)
    # Local camera covers the complete original group in BOTH actual trajectories.
    pts = np.concatenate([s["positions_m"][:, selected].reshape(-1, 3) for s in states])
    lo, hi = pts.min(0) - 2 * radius, pts.max(0) + 2 * radius
    origin = (lo + hi) / 2
    half = max(hi - lo) / 2
    bounds = np.array([-half, half] * 3) * 1e6
    wall = pv.read(wall_path())
    plot = pv.Plotter(
        off_screen=True, shape=(2, 2), row_weights=[1, 1.5], window_size=(1920, 1440), border=False
    )
    actors = []
    for col, s in enumerate(states):
        label = "Frozen Stage4L" if col == 0 else "R3"
        plot.subplot(0, col)
        decorate(
            plot, wall, (np.array(wall.bounds) * 1e6).tolist(), np.zeros(3), True, label + " | Full vessel"
        )
        global_actor = glyph(plot, s["positions_m"][0][masks[col]] * 1e6, radius * 1e6, BLUE)
        plot.subplot(1, col)
        decorate(plot, wall, bounds.tolist(), origin, True, label + " | Local view")
        local_actor = glyph(plot, (s["positions_m"][0][masks[col]] - origin) * 1e6, radius * 1e6, BLUE)
        actors.append((global_actor, local_actor))
    out = BASE / "videos"
    out.mkdir(exist_ok=True)
    stem = case + "_old_new_full_5ms_actual_radius"
    filename = out / (stem + "_60fps.mp4")
    clock = np.linspace(3, 3.005, 601)
    mapping = []
    writer = imageio.get_writer(
        filename,
        fps=60,
        codec="libx264",
        quality=8,
        pixelformat="yuv420p",
        macro_block_size=1,
        ffmpeg_params=["-movflags", "+faststart"],
    )
    try:
        for frame, t in enumerate(clock):
            row = dict(frame=frame, requested_time_s=float(t), states=[])
            for col, s in enumerate(states):
                k = max(0, min(len(s["time_s"]) - 1, int(np.searchsorted(s["time_s"], t, side="right") - 1)))
                # Only rounding-sized timestamp disagreement is admitted.
                nxt = min(k + 1, len(s["time_s"]) - 1)
                if abs(s["time_s"][nxt] - t) <= 16 * np.finfo(float).eps * max(1, abs(t)):
                    k = nxt
                x = s["positions_m"][k][masks[col]]
                actors[col][0].points = x * 1e6
                actors[col][1].points = (x - origin) * 1e6
                plot.subplot(1, col)
                label = "Frozen Stage4L" if col == 0 else "R3"
                plot.add_text(
                    f"{label} | t = {s['time_s'][k] * 1000:.6f} ms\n{len(x)} active bubbles | physical radius 1.78 um",
                    position="upper_left",
                    font_size=16,
                    color="#d5dbe3",
                    name="title",
                )
                row["states"].append(
                    dict(
                        method=label,
                        actual_time_s=float(s["time_s"][k]),
                        index=k,
                        active_count=len(x),
                        sample_age_s=float(t - s["time_s"][k]),
                    )
                )
            writer.append_data(plot.screenshot(return_img=True))
            mapping.append(row)
            if frame % 120 == 0:
                print(dict(case=case, frame=frame, of=601), flush=True)
    finally:
        writer.close()
        plot.close()
    source = stem + "_RENDER_SOURCE.py"
    (out / source).write_bytes(Path(__file__).read_bytes())
    write(
        out / (stem + "_FRAME_MAP.json"),
        dict(
            frames=mapping,
            fps=60,
            bubble_display_scale=1,
            physical_radius_m=radius,
            local_origin_m=origin.tolist(),
            local_bounds_um=bounds.tolist(),
            camera="Fixed identical old/new cameras and equal spatial scale; local camera uses union of original group extents. Every active particle remains in both scenes.",
            trajectory_interpolation=False,
            source_paths=[str(p) for p in paths],
            source_sha256=[sha(p) for p in paths],
            source_lineage="Fresh frozen Stage4L reference, separately audited; not labeled archived Stage4L",
            actual_group_IDs=group.tolist(),
        ),
    )
    audit(stem, out, source)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", choices=["B2", "B3"], required=True)
    main(p.parse_args().case)

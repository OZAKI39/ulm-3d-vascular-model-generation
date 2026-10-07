"""Uniform 60 fps video of actual accepted B2 states; real radii, no interpolation."""

import os
from pathlib import Path

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import imageio.v2 as imageio
import numpy as np
import pyvista as pv
from scripts.stage5a_r3.common import BASE, write, sha
from scripts.stage5a_r3.figures import event_data, BLUE, RED
from scripts.stage5a_r2.motion_visuals import wall_path, glyph, decorate


def trace_grid(plot, low, high):
    color = "#d5dbe3"
    span = max(high - low)
    lines = []
    for axis in [0, 1]:
        for value in np.linspace(low[axis], high[axis], 4):
            a = low.copy()
            b = low.copy()
            a[axis] = b[axis] = value
            b[1 - axis] = high[1 - axis]
            lines.extend([a, b])
    plot.add_lines(np.array(lines), color="#435060", width=0.7, connected=False)
    for axis in range(3):
        a = low.copy()
        b = a.copy()
        b[axis] = high[axis]
        plot.add_lines(np.array([a, b]), color=color, width=1.2)
        label = (a + b) / 2
        label[(axis + 1) % 3] -= 0.12 * span
        plot.add_point_labels(
            label[None],
            ["d" + "xyz"[axis] + " (nm)"],
            font_size=14,
            text_color=color,
            shape=None,
            show_points=False,
            always_visible=True,
        )
    plot.camera_position = [
        (low + high) / 2 + span * np.array([1.0, -1.55, 0.8]),
        (low + high) / 2,
        [0, 0, 1],
    ]
    plot.enable_parallel_projection()
    plot.camera.parallel_scale = 0.85 * span


def main():
    old, new, path = event_data()
    out = BASE / "videos"
    out.mkdir(parents=True, exist_ok=True)
    movie = out / "B2_actual_wall_stall_crossing_60fps.mp4"
    radius = float(new["radii_m"][np.flatnonzero(new["particle_id"] == 152)[0]])
    indices = [int(np.flatnonzero(d["particle_id"] == 152)[0]) for d in [old, new]]
    origin = new["positions_m"][0, indices[1]]
    wall = pv.read(wall_path())
    plot = pv.Plotter(
        off_screen=True, shape=(2, 3), row_weights=[1, 2.6], window_size=(1920, 1080), border=False
    )
    actors = []
    for col, data in enumerate([old, new]):
        plot.subplot(0, col)
        global_bounds = (np.array(wall.bounds) * 1e6).tolist()
        decorate(
            plot, wall, global_bounds, np.zeros(3), True, "Full vessel | " + ("R2" if col == 0 else "R3")
        )
        global_points = glyph(plot, data["positions_m"][0][data["active"][0]] * 1e6, radius * 1e6, BLUE)
        plot.subplot(1, col)
        decorate(
            plot,
            wall,
            [-5, 5, -5, 5, -5, 5],
            origin,
            True,
            "R2: original stall" if col == 0 else "R3: continued motion",
        )
        points = glyph(plot, (data["positions_m"][0][data["active"][0]] - origin) * 1e6, radius * 1e6, BLUE)
        target = glyph(
            plot, (data["positions_m"][0, indices[col]][None] - origin) * 1e6, radius * 1e6, "#75a8ce"
        )
        actors.append((global_points, points, target))
    plot.subplot(0, 2)
    plot.set_background("#080b12")
    plot.add_text(
        "Actual radius: 1.78 um\n186 active bubbles simulated\nLight blue: bubble 152\n\nStored accepted states\nNo position interpolation",
        position="upper_left",
        font_size=17,
        color="#d5dbe3",
    )
    plot.subplot(1, 2)
    plot.set_background("#080b12")
    trace = [(d["positions_m"][:, i] - origin) * 1e9 for d, i in zip([old, new], indices)]
    combined = np.concatenate(trace)
    mid = (combined.max(0) + combined.min(0)) / 2
    span = np.ptp(combined, axis=0).max() * 1.1
    low, high = mid - span / 2, mid + span / 2
    trace_grid(plot, low, high)
    plot.add_text("Center displacement (nm)", position="upper_left", font_size=16, color="#d5dbe3")
    for x, color in zip(trace, [RED, BLUE]):
        plot.add_mesh(pv.lines_from_points(x), color=color, line_width=2, opacity=0.45)
    markers = [glyph(plot, x[0:1], span * 0.016, color) for x, color in zip(trace, [RED, BLUE])]
    clock = np.linspace(float(new["time_s"][0]), float(new["time_s"][-1]), 361)
    mapping = []
    writer = imageio.get_writer(
        movie,
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
            for col, data in enumerate([old, new]):
                k = max(
                    0, min(len(data["time_s"]) - 1, int(np.searchsorted(data["time_s"], t, side="right") - 1))
                )
                x = data["positions_m"][k]
                mask = data["active"][k]
                assert np.array_equal(mask, data["active"][0])
                actors[col][0].points = x[mask] * 1e6
                actors[col][1].points = (x[mask] - origin) * 1e6
                actors[col][2].points = (x[indices[col]][None] - origin) * 1e6
                markers[col].points = trace[col][k : k + 1]
                plot.subplot(1, col)
                stalled = col == 0 and t > old["time_s"][-1]
                label = "R2: stopped" if stalled else "R2" if col == 0 else "R3"
                plot.add_text(
                    f"{label} | t = {data['time_s'][k] * 1000:.6f} ms\nLocal view | actual spheres",
                    position="upper_left",
                    font_size=15,
                    color="#d5dbe3",
                    name="title",
                )
                row["states"].append(
                    dict(
                        method="R2" if col == 0 else "R3",
                        index=k,
                        actual_time_s=float(data["time_s"][k]),
                        active_count=int(mask.sum()),
                        stopped=stalled,
                    )
                )
            writer.append_data(plot.screenshot(return_img=True))
            mapping.append(row)
    finally:
        writer.close()
        plot.close()
    write(
        out / "B2_actual_wall_stall_crossing_FRAME_MAP.json",
        dict(
            frames=mapping,
            fps=60,
            bubble_display_scale=1.0,
            physical_radius_m=radius,
            local_origin_m=origin.tolist(),
            simulation_source=str(path),
            simulation_sha256=sha(path),
            scalar_trace_scope="Actual bubble-center displacement with separate labelled nm axes; markers are plotting symbols, not rescaled bubbles",
            trajectory_interpolation=False,
            old_stop_policy="Hold last recorded R2 state and explicitly display stopped status",
        ),
    )
    from scripts.stage5a_r3.video_audit import audit

    (out / "B2_RENDER_SOURCE.py").write_bytes(Path(__file__).read_bytes())
    audit()
    print(movie, flush=True)


if __name__ == "__main__":
    main()

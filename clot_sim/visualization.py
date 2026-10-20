"""Fixed-camera replay of real saved states, with red phases and world Z upright.

Coordinates and surviving bonds are shown exactly. Missing intermediate states
are never interpolated; a sparse input selection is labelled as a sparse replay.
"""

from pathlib import Path
from itertools import product
import argparse
import json
import subprocess
import shutil
import imageio_ffmpeg
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from matplotlib.ticker import MaxNLocator
from .config import INPUTS, DX
from .io import sha, dump, load_npz, verify_inputs
from .geometry import load_reference
from .fragmentation import components

FFMPEG = shutil.which("ffmpeg") or imageio_ffmpeg.get_ffmpeg_exe()

W, H, FPS = 2560, 2048, 30
VIEW = (180, 140, 2180, 1640)
AZIMUTH, ELEVATION = -68.0, -16.0
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONTS = {s: ImageFont.truetype(FONT, s) for s in (36, 40, 46, 52, 60)}
C = dict(
    weak="#ca5360",
    strong="#8b283e",
    anchor="#8fa4b0",
    bond="#865662",
    interface="#8198a3",
    grid="#8c9aa4",
    text="#35444d",
)


def basis():
    az, el = np.deg2rad([AZIMUTH, ELEVATION])
    out = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    right = np.cross(-out, [0.0, 0.0, 1.0])
    right /= np.linalg.norm(right)
    up = np.cross(right, -out)
    return out, right, up


def lines(points, edges):
    used, inverse = np.unique(edges, return_inverse=True)
    mesh = pv.PolyData(points[used])
    mesh.verts = np.empty(0, dtype=np.int64)
    mesh.lines = np.c_[np.full(len(edges), 2), inverse.reshape(-1, 2)].ravel()
    return mesh


class Scene:
    def __init__(self, bounds, kind):
        self.bounds, self.kind = bounds, kind
        self.focal = bounds.mean(1)
        self.out, self.right, self.up = basis()
        self.corners = np.array(list(product(*bounds)))
        q = self.corners - self.focal
        self.scale = (
            max(np.ptp(q @ self.up) / 2, np.ptp(q @ self.right) / 2 * VIEW[3] / VIEW[2])
            * 1.04
        )
        self.ticks = []
        for lo, hi in bounds:
            values = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
            self.ticks.append(values[(values > lo + 1e-9) & (values < hi - 1e-9)])
        if kind == "closeup":
            self.ticks = [
                np.array([-6.6, -6.0, -5.4]),
                np.array([-0.5, 0.0, 0.5]),
                np.array([0.6, 1.0, 1.4]),
            ]
        self.p = pv.Plotter(
            off_screen=True, window_size=VIEW[2:], lighting="three lights"
        )
        self.p.set_background("white")
        self.p.render_window.SetMultiSamples(0)
        self.p.render_window.SetAlphaBitPlanes(1)
        self.p.enable_depth_peeling(number_of_peels=8, occlusion_ratio=0.01)
        grid = []
        for axis in range(3):
            other = [j for j in range(3) if j != axis]
            for side in bounds[axis]:
                for along, across in [other, other[::-1]]:
                    for value in np.unique(np.r_[bounds[across], self.ticks[across]]):
                        edge = np.zeros((2, 3))
                        edge[:, axis], edge[:, across], edge[:, along] = (
                            side,
                            value,
                            bounds[along],
                        )
                        grid.append(edge)
        grid = np.array(list({tuple(e.ravel()): e for e in grid}.values()))
        self.grid_faces = [
            int(np.isclose(grid[:, :, axis], side).all(1).sum())
            for axis in range(3)
            for side in bounds[axis]
        ]
        assert min(self.grid_faces) >= 8
        self.p.add_mesh(
            lines(grid.reshape(-1, 3), np.arange(2 * len(grid)).reshape(-1, 2)),
            color=C["grid"],
            opacity=0.70,
            line_width=1.6,
            lighting=False,
        )
        self.p.enable_parallel_projection()
        self.p.camera_position = [self.focal + 50 * self.out, self.focal, self.up]
        self.p.camera.parallel_scale = self.scale
        self.p.camera.clipping_range = (0.1, 100.0)
        self.fixed_camera = np.array(list(self.p.camera_position), dtype=float)
        self.axis_edges = self.edges()
        self.axis_overlay, self.label_boxes = self.axes()

    def project(self, point):
        q = point - self.focal
        return np.array(
            [
                VIEW[0] + VIEW[2] / 2 + (q @ self.right) * VIEW[3] / (2 * self.scale),
                VIEW[1] + VIEW[3] / 2 - (q @ self.up) * VIEW[3] / (2 * self.scale),
            ]
        )

    def edges(self):
        xy = np.c_[self.corners @ self.right, self.corners @ self.up]
        center = xy.mean(0)
        result = []
        for axis in range(3):
            other = [j for j in range(3) if j != axis]
            candidates = []
            for vals in product(*[self.bounds[j] for j in other]):
                edge = np.zeros((2, 3))
                edge[:, axis] = self.bounds[axis]
                for j, value in zip(other, vals):
                    edge[:, j] = value
                q = np.c_[edge @ self.right, edge @ self.up]
                d = q[1] - q[0]
                sides = d[0] * (xy[:, 1] - q[0, 1]) - d[1] * (xy[:, 0] - q[0, 0])
                if not (sides.min() >= -1e-10 or sides.max() <= 1e-10):
                    continue
                normal = np.array([-d[1], d[0]])
                normal /= np.linalg.norm(normal)
                if normal @ (q.mean(0) - center) < 0:
                    normal = -normal
                score = q.mean(0)[0] if axis == 2 else -q.mean(0)[1]
                candidates.append((score, edge, normal * [1.0, -1.0]))
            _, edge, normal = max(candidates, key=lambda a: a[0])
            result.append((axis, edge, normal))
        return result

    def axes(self):
        layer = Image.new("RGBA", (W, H))
        draw = ImageDraw.Draw(layer)
        boxes = []
        for axis, edge, normal in self.axis_edges:
            xy = np.array([self.project(p) for p in edge])
            draw.line([tuple(q) for q in xy], fill=C["grid"], width=2)
            for value in self.ticks[axis]:
                p = edge.mean(0)
                p[axis] = value
                xyv = self.project(p)
                draw.line(
                    [tuple(xyv), tuple(xyv + normal * 11)], fill=C["text"], width=2
                )
                pos = xyv + normal * (86 if axis == 2 else 58)
                text = f"{value:g}"
                draw.text(tuple(pos), text, font=FONTS[52], fill=C["text"], anchor="mm")
                boxes.append(
                    (text, draw.textbbox(tuple(pos), text, font=FONTS[52], anchor="mm"))
                )
            pos = (
                xy[np.argmin(xy[:, 1])] + [0.0, -64.0]
                if axis == 2
                else self.project(edge.mean(0)) + normal * 158
            )
            text = "XYZ"[axis] + " (mm)"
            draw.text(tuple(pos), text, font=FONTS[60], fill=C["text"], anchor="mm")
            boxes.append(
                (text, draw.textbbox(tuple(pos), text, font=FONTS[60], anchor="mm"))
            )
        for text, b in boxes:
            assert b[0] > 0 and b[1] > 100 and b[2] < W and b[3] < H - 100, (
                self.kind,
                text,
                b,
            )
        for i, (a, ba) in enumerate(boxes):
            for b, bb in boxes[i + 1 :]:
                assert not (
                    min(ba[2], bb[2]) - max(ba[0], bb[0]) > 2
                    and min(ba[3], bb[3]) - max(ba[1], bb[1]) > 2
                ), (a, b)
        z0, z1 = self.project(self.focal), self.project(self.focal + [0.0, 0.0, 0.2])
        assert abs(z1[0] - z0[0]) < 1e-8 and z1[1] < z0[1]
        return layer, boxes

    def update(self, ref, field, s):
        x, n = s["x"] * 1e3, int(ref["clot_count"])
        visible = ((x >= self.bounds[:, 0]) & (x <= self.bounds[:, 1])).all(1)
        if self.kind == "full":
            assert visible.all()
        bonds, typ = ref["bonds"], ref["bond_types"]
        alive = s["d"] < 1.0
        counts = {}
        for t, color in [(1, C["bond"]), (2, C["interface"])]:
            selected = alive & (typ == t) & visible[bonds].all(1)
            # Use true reference distances, including jittered short neighbours.
            short = selected & (ref["r0"] <= 1.35 * DX)
            for label, mask, alpha, width in [
                ("all", selected, 0.045, 0.8),
                ("short", short, 0.25 if t == 1 else 0.16, 1.2),
            ]:
                name = f"bonds_{t}_{label}"
                if mask.any():
                    self.p.add_mesh(
                        lines(x, bonds[mask]),
                        color=color,
                        opacity=alpha,
                        line_width=width,
                        lighting=False,
                        name=name,
                        reset_camera=False,
                    )
                else:
                    self.p.remove_actor(name, reset_camera=False)
            counts[str(t)] = dict(
                visible_alive=int(selected.sum()), emphasized=int(short.sum())
            )
            assert not np.any(selected & ~alive)
        for name, ids, color, size, alpha in [
            ("anchors", np.arange(n, len(x)), C["anchor"], 7.0, 0.55),
            (
                "weak",
                np.flatnonzero(field["material_phase"][:n] == 0),
                C["weak"],
                11.0,
                1.0,
            ),
            (
                "strong",
                np.flatnonzero(field["material_phase"][:n] == 1),
                C["strong"],
                11.0,
                1.0,
            ),
        ]:
            ids = ids[visible[ids]]
            if len(ids):
                self.p.add_points(
                    x[ids],
                    color=color,
                    point_size=size,
                    opacity=alpha,
                    render_points_as_spheres=True,
                    ambient=0.45,
                    diffuse=0.55,
                    specular=0.10,
                    name=name,
                    reset_camera=False,
                )
            else:
                self.p.remove_actor(name, reset_camera=False)
        assert np.allclose(
            np.array(list(self.p.camera_position), dtype=float),
            self.fixed_camera,
            rtol=0,
            atol=1e-12,
        )
        self.p.render()
        image = Image.new("RGBA", (W, H), "white")
        image.paste(
            Image.fromarray(self.p.screenshot(return_img=True)).convert("RGBA"),
            VIEW[:2],
        )
        image.alpha_composite(self.axis_overlay)
        draw = ImageDraw.Draw(image)
        draw.text(
            (W - 130, 58),
            f't = {float(s["time"])*1e3:.5f} ms',
            font=FONTS[46],
            fill=C["text"],
            anchor="rt",
        )
        outside = int((~visible[:n]).sum())
        if outside:
            draw.text(
                (130, 58),
                f"Outside close-up: {outside} clot nodes",
                font=FONTS[40],
                fill=C["text"],
                anchor="lt",
            )
        items = [
            ("dot", "Weak phase", C["weak"]),
            ("dot", "Strong phase", C["strong"]),
            ("dot", "Fixed anchors", C["anchor"]),
            ("line", "Surviving bonds", C["bond"]),
        ]
        widths = [48 + draw.textlength(label, font=FONTS[40]) for _, label, _ in items]
        gap = 65
        xx, yy = (W - sum(widths) - gap * (len(items) - 1)) / 2, H - 55
        for (kind, label, color), width in zip(items, widths):
            if kind == "dot":
                draw.ellipse((xx + 3, yy - 9, xx + 21, yy + 9), fill=color)
            else:
                draw.line([(xx - 3, yy), (xx + 26, yy)], fill=color, width=3)
            draw.text((xx + 48, yy), label, font=FONTS[40], fill=C["text"], anchor="lm")
            xx += width + gap
        return image.convert("RGB"), dict(
            outside_clot_nodes=outside,
            visible_clot_nodes=int(visible[:n].sum()),
            bonds=counts,
            camera_unchanged=True,
            broken_edges_drawn=0,
        )


def encoder(path):
    assert not path.exists(), f"Refuse to overwrite finalized video: {path}"
    return subprocess.Popen(
        [
            FFMPEG,
            "-v",
            "error",
            "-f",
            "rawvideo",
            "-pix_fmt",
            "rgb24",
            "-s",
            f"{W}x{H}",
            "-r",
            str(FPS),
            "-i",
            "-",
            "-an",
            "-c:v",
            "libx264",
            "-crf",
            "17",
            "-preset",
            "medium",
            "-threads",
            "4",
            "-pix_fmt",
            "yuv420p",
            "-movflags",
            "+faststart",
            str(path),
        ],
        stdin=subprocess.PIPE,
    )


def validate_video(path, frames):
    reader = imageio_ffmpeg.read_frames(str(path), pix_fmt="rgb24")
    metadata = next(reader)
    reader.close()
    decoded = subprocess.run(
        [
            FFMPEG,
            "-v",
            "error",
            "-xerror",
            "-i",
            str(path),
            "-progress",
            "pipe:1",
            "-f",
            "null",
            "-",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    counts = [
        int(line.split("=", 1)[1])
        for line in decoded.stdout.splitlines()
        if line.startswith("frame=")
    ]
    assert counts and counts[-1] == frames, (counts, frames)
    return dict(
        file=str(path),
        sha256=sha(path),
        width=metadata["size"][0],
        height=metadata["size"][1],
        frames=frames,
        fps=metadata["fps"],
        duration_s=metadata["duration"],
        full_decode="PASS",
    )


def render(paths, output, video=False):
    verify_inputs()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    ref = load_reference()
    field = load_npz(INPUTS / "material.npz")
    unique = {}
    for path in paths:
        path = Path(path)
        with np.load(path, allow_pickle=False) as z:
            if all(k in z for k in ["x", "d", "time", "step", "atom_ids"]):
                unique[float(z["time"])] = path
    if not unique:
        raise ValueError("No saved x/d/time/step/atom_ids states found")
    ordered = sorted(unique.items())
    states = []
    sources = []
    previous = np.zeros(len(ref["bonds"]))
    for time, path in ordered:
        z = load_npz(path)
        state = {k: z[k] for k in ["x", "d", "time", "step", "atom_ids"]}
        assert np.array_equal(state["atom_ids"], ref["atom_ids"])
        assert np.all(state["d"] >= previous - 1e-12)
        assert np.isfinite(state["x"]).all()
        assert np.array_equal(
            state["x"][int(ref["clot_count"]) :], ref["X"][int(ref["clot_count"]) :]
        )
        previous = state["d"]
        states.append(state)
        sources.append(dict(path=str(path.resolve()), sha256=sha(path)))
    xyz = np.concatenate([s["x"] * 1e3 for s in states])
    bounds = np.column_stack((xyz.min(0) - 0.12, xyz.max(0) + 0.12))
    scene = Scene(bounds, "full")
    sink = None
    records = []
    try:
        if video:
            sink = encoder(output / "replay.mp4")
        for i, (state, source) in enumerate(zip(states, sources)):
            image, detail = scene.update(ref, field, state)
            image.save(output / f"state_{i:04d}.png")
            top = components(ref, state["d"])
            records.append(
                dict(
                    **source,
                    time_s=float(state["time"]),
                    step=int(state["step"]),
                    free_components=int((~top["attached"]).sum()),
                    **detail,
                )
            )
            if sink:
                for _ in range(FPS):
                    sink.stdin.write(np.asarray(image).tobytes())
        if sink:
            sink.stdin.close()
            assert sink.wait() == 0
        audit = dict(
            saved_states=len(states),
            time_range_s=[ordered[0][0], ordered[-1][0]],
            interpolated=False,
            timing="One second per supplied state; timestamps show actual simulated time",
            world_Z_vertical=True,
            camera_rotation=False,
            grid_segments_per_face=scene.grid_faces,
            frames=records,
        )
        if video:
            audit["video"] = validate_video(output / "replay.mp4", FPS * len(states))
        dump(output / "RENDER_AUDIT.json", audit)
    finally:
        scene.p.close()
        if sink and sink.poll() is None:
            sink.kill()
    for source in sources:
        assert sha(source["path"]) == source["sha256"]
    (output / "index.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Clot saved-state replay</title><h1>Saved-state replay</h1><p>Real coordinates and bonds; one second per supplied state; no interpolation.</p>'
        + (' <video controls width="900" src="replay.mp4"></video>' if video else "")
        + "".join(
            f'<p><img width="900" src="state_{i:04d}.png"></p>'
            for i in range(len(states))
        )
    )
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--video", action="store_true")
    args = parser.parse_args()
    paths = []
    for path in args.states:
        paths.extend(sorted(path.glob("*.npz")) if path.is_dir() else [path])
    render(paths, args.output, args.video)


if __name__ == "__main__":
    main()

"""Fixed-camera replay of accepted V4 fracture states, without a physics rerun."""
from pathlib import Path
from itertools import product
import argparse
import hashlib
import json
import subprocess
import time

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont, ImageColor
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from matplotlib.ticker import MaxNLocator

R = Path(__file__).resolve().parents[3]
TASK = "y_clot_fracture_reference_style_v2"
RUN = "20261020T094240Z"
SOURCE = R / "verification/y_clot_irregular_current_method_v4/20261020T041100Z"
CASE = SOURCE / "cases/H_Y_IRREGULAR_LPS_V4"
OUT = R / "results" / TASK / RUN
AUDIT = R / "verification" / TASK / RUN
W, H, FPS = 2560, 2048, 30
VIEW = (180, 140, 2180, 1640)
AZIMUTH, ELEVATION = -68.0, -16.0
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONTS = {s: ImageFont.truetype(FONT, s) for s in (36, 40, 46, 52, 60)}
C = dict(weak="#ca5360", strong="#8b283e", anchor="#8fa4b0",
         bond="#865662", interface="#8198a3", grid="#8c9aa4", text="#35444d")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(path, value):
    def convert(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, Path):
            return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, default=convert) + "\n")


def load(path, keys=None):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in (z.files if keys is None else keys)}


def inputs():
    ref = load(SOURCE / "provenance/reference.npz")
    field = load(SOURCE / "BIPHASIC_FIELD.npz")
    paths = [CASE / "keys/initial.npz", *sorted((CASE / "keys").glob("frame_*.npz"))]
    for name in ("pilot", "final", "first_softening", "first_failure",
                 "first_free_component", "first_target_observed", "max_target_count"):
        p = CASE / "keys" / (name + ".npz")
        if p.exists():
            paths.append(p)
    if (CASE / "last_accepted.npz").exists():
        paths.append(CASE / "last_accepted.npz")
    unique = {}
    for p in paths:
        with np.load(p) as z:
            unique[float(z["time"])] = p
    paths = sorted(unique.items())
    result = json.loads((SOURCE / "NATIVE_V4_RESULT.json").read_text())
    assert result["status"] == "COMPLETED"
    assert abs(paths[-1][0] - result["accepted_time_s"]) < 1e-14
    states = [load(p, ("x", "d", "atom_ids", "time", "step")) for _, p in paths]
    n = int(ref["clot_count"])
    assert n == 1927 and len(ref["X"]) == 3488 and len(ref["bonds"]) == 153986
    assert np.array_equal(states[0]["x"], ref["X"])
    last_d = states[0]["d"]
    for s in states:
        assert np.array_equal(s["atom_ids"], field["atom_ids"])
        assert np.array_equal(s["x"][n:], ref["X"][n:])
        assert np.all(s["d"] >= last_d - 1e-12)
        assert np.isfinite(s["x"]).all()
        last_d = s["d"]
    sources = [SOURCE / "provenance/reference.npz", SOURCE / "BIPHASIC_FIELD.npz",
               SOURCE / "NATIVE_V4_RESULT.json", SOURCE / "IRREGULAR_REFERENCE_GEOMETRY.json"]
    sources += [p for _, p in paths]
    hashes = {str(p.relative_to(R)): sha(p) for p in sources}
    frozen = AUDIT / "SOURCE_HASHES.json"
    if frozen.exists():
        assert json.loads(frozen.read_text()) == hashes
    else:
        dump(frozen, hashes)
    return ref, field, paths, states, hashes


def basis():
    az, el = np.deg2rad([AZIMUTH, ELEVATION])
    out = np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])
    right = np.cross(-out, [0., 0., 1.])
    right /= np.linalg.norm(right)
    up = np.cross(right, -out)
    return out, right, up


def lines(points, edges):
    used, inverse = np.unique(edges, return_inverse=True)
    mesh = pv.PolyData(points[used])
    mesh.verts = np.empty(0, dtype=np.int64)
    mesh.lines = np.c_[np.full(len(edges), 2), inverse.reshape(-1, 2)].ravel()
    return mesh


def graph_stats(ref, s):
    n, N = int(ref["clot_count"]), len(ref["X"])
    b = ref["bonds"][s["d"] < 1]
    nc, labels = connected_components(
        csr_matrix((np.ones(2*len(b)),
                    (np.r_[b[:, 0], b[:, 1]], np.r_[b[:, 1], b[:, 0]])), shape=(N, N)),
        directed=False)
    attached = np.zeros(nc, bool)
    attached[labels[n:]] = True
    sizes = np.bincount(labels)
    return dict(free_nodes=int((~attached[labels[:n]]).sum()),
                free_components=int((~attached).sum()),
                target_20_200=int((~attached & (sizes >= 20) & (sizes <= 200)).sum()),
                failed_bulk=int(((s["d"] >= 1) & (ref["bond_types"] == 1)).sum()),
                failed_interface=int(((s["d"] >= 1) & (ref["bond_types"] == 2)).sum()))


class Scene:
    def __init__(self, bounds, kind):
        self.bounds, self.kind = bounds, kind
        self.focal = bounds.mean(1)
        self.out, self.right, self.up = basis()
        self.corners = np.array(list(product(*bounds)))
        q = self.corners - self.focal
        self.scale = max(np.ptp(q @ self.up) / 2,
                         np.ptp(q @ self.right) / 2 * VIEW[3] / VIEW[2]) * 1.04
        self.ticks = []
        for lo, hi in bounds:
            values = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
            self.ticks.append(values[(values > lo+1e-9) & (values < hi-1e-9)])
        if kind == "closeup":
            self.ticks = [np.array([-6.6, -6., -5.4]), np.array([-.5, 0., .5]),
                          np.array([.6, 1., 1.4])]
        self.p = pv.Plotter(off_screen=True, window_size=VIEW[2:], lighting="three lights")
        self.p.set_background("white")
        self.p.render_window.SetMultiSamples(0)
        self.p.render_window.SetAlphaBitPlanes(1)
        self.p.enable_depth_peeling(number_of_peels=8, occlusion_ratio=.01)
        grid = []
        for axis in range(3):
            other = [j for j in range(3) if j != axis]
            for side in bounds[axis]:
                for along, across in [other, other[::-1]]:
                    for value in np.unique(np.r_[bounds[across], self.ticks[across]]):
                        edge = np.zeros((2, 3))
                        edge[:, axis], edge[:, across], edge[:, along] = side, value, bounds[along]
                        grid.append(edge)
        grid = np.array(list({tuple(e.ravel()): e for e in grid}.values()))
        self.grid_faces = [int(np.isclose(grid[:, :, axis], side).all(1).sum())
                           for axis in range(3) for side in bounds[axis]]
        assert min(self.grid_faces) >= 8
        self.p.add_mesh(lines(grid.reshape(-1, 3), np.arange(2*len(grid)).reshape(-1, 2)),
                        color=C["grid"], opacity=.70, line_width=1.6, lighting=False)
        self.p.enable_parallel_projection()
        self.p.camera_position = [self.focal + 50*self.out, self.focal, self.up]
        self.p.camera.parallel_scale = self.scale
        self.p.camera.clipping_range = (.1, 100.)
        self.fixed_camera = np.array(list(self.p.camera_position), dtype=float)
        self.axis_edges = self.edges()
        self.axis_overlay, self.label_boxes = self.axes()

    def project(self, point):
        q = point-self.focal
        return np.array([VIEW[0]+VIEW[2]/2+(q @ self.right)*VIEW[3]/(2*self.scale),
                         VIEW[1]+VIEW[3]/2-(q @ self.up)*VIEW[3]/(2*self.scale)])

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
                d = q[1]-q[0]
                sides = d[0]*(xy[:, 1]-q[0, 1])-d[1]*(xy[:, 0]-q[0, 0])
                if not (sides.min() >= -1e-10 or sides.max() <= 1e-10):
                    continue
                normal = np.array([-d[1], d[0]])
                normal /= np.linalg.norm(normal)
                if normal @ (q.mean(0)-center) < 0:
                    normal = -normal
                score = q.mean(0)[0] if axis == 2 else -q.mean(0)[1]
                candidates.append((score, edge, normal * [1., -1.]))
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
                draw.line([tuple(xyv), tuple(xyv+normal*11)], fill=C["text"], width=2)
                pos = xyv+normal*(86 if axis == 2 else 58)
                text = f"{value:g}"
                draw.text(tuple(pos), text, font=FONTS[52], fill=C["text"], anchor="mm")
                boxes.append((text, draw.textbbox(tuple(pos), text, font=FONTS[52], anchor="mm")))
            pos = (xy[np.argmin(xy[:, 1])] + [0., -64.] if axis == 2 else
                   self.project(edge.mean(0)) + normal*158)
            text = "XYZ"[axis] + " (mm)"
            draw.text(tuple(pos), text, font=FONTS[60], fill=C["text"], anchor="mm")
            boxes.append((text, draw.textbbox(tuple(pos), text, font=FONTS[60], anchor="mm")))
        for text, b in boxes:
            assert b[0] > 0 and b[1] > 100 and b[2] < W and b[3] < H-100, (self.kind, text, b)
        for i, (a, ba) in enumerate(boxes):
            for b, bb in boxes[i+1:]:
                assert not (min(ba[2], bb[2])-max(ba[0], bb[0]) > 2 and
                            min(ba[3], bb[3])-max(ba[1], bb[1]) > 2), (a, b)
        z0, z1 = self.project(self.focal), self.project(self.focal+[0., 0., .2])
        assert abs(z1[0]-z0[0]) < 1e-8 and z1[1] < z0[1]
        return layer, boxes

    def update(self, ref, field, s):
        x, n = s["x"]*1e3, int(ref["clot_count"])
        visible = ((x >= self.bounds[:, 0]) & (x <= self.bounds[:, 1])).all(1)
        if self.kind == "full":
            assert visible.all()
        bonds, typ = ref["bonds"], ref["bond_types"]
        alive = s["d"] < 1.
        counts = {}
        for t, color in [(1, C["bond"]), (2, C["interface"])]:
            selected = alive & (typ == t) & visible[bonds].all(1)
            # Use true reference distances, including jittered short neighbours.
            short = selected & (ref["r0"] <= 1.35*50e-6)
            for label, mask, alpha, width in [("all", selected, .045, .8),
                                             ("short", short, .25 if t == 1 else .16, 1.2)]:
                name = f"bonds_{t}_{label}"
                if mask.any():
                    self.p.add_mesh(lines(x, bonds[mask]), color=color, opacity=alpha,
                                    line_width=width, lighting=False, name=name, reset_camera=False)
                else:
                    self.p.remove_actor(name, reset_camera=False)
            counts[str(t)] = dict(visible_alive=int(selected.sum()), emphasized=int(short.sum()))
            assert not np.any(selected & ~alive)
        for name, ids, color, size, alpha in [
            ("anchors", np.arange(n, len(x)), C["anchor"], 7., .55),
            ("weak", np.flatnonzero(field["material_phase"][:n] == 0), C["weak"], 11., 1.),
            ("strong", np.flatnonzero(field["material_phase"][:n] == 1), C["strong"], 11., 1.)]:
            ids = ids[visible[ids]]
            if len(ids):
                self.p.add_points(x[ids], color=color, point_size=size, opacity=alpha,
                                  render_points_as_spheres=True, ambient=.45, diffuse=.55,
                                  specular=.10, name=name, reset_camera=False)
            else:
                self.p.remove_actor(name, reset_camera=False)
        assert np.allclose(np.array(list(self.p.camera_position), dtype=float),
                           self.fixed_camera, rtol=0, atol=1e-12)
        self.p.render()
        image = Image.new("RGBA", (W, H), "white")
        image.paste(Image.fromarray(self.p.screenshot(return_img=True)).convert("RGBA"), VIEW[:2])
        image.alpha_composite(self.axis_overlay)
        draw = ImageDraw.Draw(image)
        draw.text((W-130, 58), f't = {float(s["time"])*1e3:.5f} ms',
                  font=FONTS[46], fill=C["text"], anchor="rt")
        outside = int((~visible[:n]).sum())
        if outside:
            draw.text((130, 58), f"Outside close-up: {outside} clot nodes",
                      font=FONTS[40], fill=C["text"], anchor="lt")
        items = [("dot", "Weak phase", C["weak"]), ("dot", "Strong phase", C["strong"]),
                 ("dot", "Fixed anchors", C["anchor"]), ("line", "Surviving bonds", C["bond"])]
        widths = [48 + draw.textlength(label, font=FONTS[40]) for _, label, _ in items]
        gap = 65
        xx, yy = (W-sum(widths)-gap*(len(items)-1))/2, H-55
        for (kind, label, color), width in zip(items, widths):
            if kind == "dot":
                draw.ellipse((xx+3, yy-9, xx+21, yy+9), fill=color)
            else:
                draw.line([(xx-3, yy), (xx+26, yy)], fill=color, width=3)
            draw.text((xx+48, yy), label, font=FONTS[40], fill=C["text"], anchor="lm")
            xx += width+gap
        return image.convert("RGB"), dict(outside_clot_nodes=outside,
                    visible_clot_nodes=int(visible[:n].sum()), bonds=counts,
                    camera_unchanged=True, broken_edges_drawn=0)


def encoder(path):
    assert not path.exists(), f"Refuse to overwrite finalized video: {path}"
    return subprocess.Popen(["ffmpeg", "-v", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
        "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-", "-an", "-c:v", "libx264",
        "-crf", "17", "-preset", "medium", "-threads", "4", "-pix_fmt", "yuv420p",
        "-movflags", "+faststart", str(path)], stdin=subprocess.PIPE)


def validate_video(path, frames):
    metadata = json.loads(subprocess.check_output(["ffprobe", "-v", "error", "-show_streams",
                                                   "-show_format", "-of", "json", str(path)]))
    v = next(s for s in metadata["streams"] if s["codec_type"] == "video")
    assert int(v["nb_frames"]) == frames and v["pix_fmt"] == "yuv420p"
    subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-i", str(path),
                    "-f", "null", "-"], check=True, stdout=subprocess.DEVNULL)
    return dict(file=str(path.relative_to(OUT)), sha256=sha(path),
                width=v["width"], height=v["height"], frames=frames,
                fps=v["r_frame_rate"], duration_s=float(v["duration"]), full_decode="PASS")


def main(preview=False):
    start = time.monotonic()
    for p in [AUDIT, OUT/"videos", OUT/"figures", OUT/"visual_review"]:
        p.mkdir(parents=True, exist_ok=True)
    ref, field, paths, states, hashes = inputs()
    all_x = np.concatenate([s["x"]*1e3 for s in states])
    full = np.column_stack((all_x.min(0)-.12, all_x.max(0)+.12))
    close = np.array([[-6.8, -5.2], [-.70, .70], [.50, 1.55]])
    assert ((states[0]["x"]*1e3 >= close[:, 0]) & (states[0]["x"]*1e3 <= close[:, 1])).all()
    scenes = {k: Scene(b, k) for k, b in [("closeup", close), ("full", full)]}
    dump(AUDIT/"CAMERAS.json", {k: dict(bounds_mm=s.bounds, position=s.fixed_camera,
        parallel_scale_mm=s.scale, grid_segments_per_face=s.grid_faces, labels=s.label_boxes,
        camera_motion=False, world_Z_vertical=True) for k, s in scenes.items()})
    requested = [0., .00054, .00112, .00185, .0024, .0036]
    keyframes = sorted(set(min(range(len(paths)), key=lambda i: abs(paths[i][0]-t)) for t in requested))
    indices = keyframes if preview else range(len(states))
    sinks, videos = {}, {}
    records, thumbnails = [], {k: [] for k in scenes}
    try:
        if not preview:
            for kind in scenes:
                path = OUT/"videos"/f"VIDEO_V4_FRACTURE_{kind.upper()}_2560.mp4"
                videos[kind], sinks[kind] = path, encoder(path)
        for i in indices:
            s, source = states[i], paths[i][1]
            stats = graph_stats(ref, s)
            record = dict(source=str(source.relative_to(R)), source_sha256=hashes[str(source.relative_to(R))],
                          time_s=float(s["time"]), step=int(s["step"]), repeat=30 if i in [0, len(states)-1] else 3,
                          **stats, views={})
            for kind, scene in scenes.items():
                image, detail = scene.update(ref, field, s)
                record["views"][kind] = detail
                if not preview:
                    raw = np.asarray(image).tobytes()
                    for _ in range(record["repeat"]):
                        sinks[kind].stdin.write(raw)
                if i in keyframes:
                    dest = OUT/("visual_review" if preview else "figures")
                    image.save(dest/f"{kind.upper()}_{i:03d}.png")
                    small = image.copy()
                    small.thumbnail((640, 512))
                    thumbnails[kind].append(small)
            records.append(record)
            if i % 10 == 0 or preview:
                print("STATE", i, "/", len(states), "time_ms", float(s["time"])*1e3,
                      "elapsed_s", round(time.monotonic()-start, 1), flush=True)
        for enc in sinks.values():
            enc.stdin.close()
        for enc in sinks.values():
            assert enc.wait() == 0
    finally:
        for scene in scenes.values():
            scene.p.close()
        for enc in sinks.values():
            if enc.poll() is None:
                enc.kill()
    for kind, images in thumbnails.items():
        sheet = Image.new("RGB", (1920, 1024), "white")
        for i, image in enumerate(images):
            sheet.paste(image, ((i % 3)*640, (i // 3)*512))
        sheet.save(OUT/"visual_review"/f"{kind.upper()}_CONTACT_SHEET.jpg", quality=95)
    for p, digest in hashes.items():
        assert sha(R/p) == digest
    dump(AUDIT/("PREVIEW_FRAMES.json" if preview else "VIDEO_FRAMES.json"), records)
    if preview:
        dump(AUDIT/"PREVIEW.json", dict(status="PASS", real_saved_states=True,
                                      physical_time_range_s=[paths[0][0], paths[-1][0]]))
        print("PREVIEW_COMPLETE", flush=True)
        return
    frames = sum(q["repeat"] for q in records)
    checks = []
    for path in videos.values():
        checks.append(validate_video(path, frames))
        ppt = path.with_name(path.name.replace("_2560", "_PPT_1600"))
        assert not ppt.exists()
        subprocess.run(["ffmpeg", "-v", "error", "-i", str(path), "-vf", "scale=1600:1280:flags=lanczos",
                        "-c:v", "libx264", "-crf", "17", "-preset", "medium", "-threads", "4",
                        "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(ppt)], check=True)
        checks.append(validate_video(ppt, frames))
    dump(AUDIT/"VALIDATION.json", dict(status="PASS", source_hashes_unchanged=True,
        states=len(states), frames=frames, physical_time_s=[paths[0][0], paths[-1][0]],
        world_Z_vertical=True, camera_rotation=False, camera_translation=False,
        physical_displacement_scale=1, position_interpolation=False, invented_nodes=0,
        broken_edges_drawn=0, full_view_missing_nodes=0, current_material_phases_preserved=True,
        video_timing="Each saved state held for 3 frames; first/last held for 30. Actual state timestamp shown.",
        scientific_scope="Accepted trajectory replay, not a new physics run. Markers indicate positions/phases, not radii.",
        outputs=checks, renderer_sha256=sha(__file__), wall_s=time.monotonic()-start))
    html = """<!doctype html><meta charset="utf-8"><title>V4 血栓真实碎裂过程</title>
    <style>body{font:18px system-ui;max-width:1100px;margin:35px auto;color:#35444d}
    video,img{width:100%;background:white}a{color:#9a2741}p{line-height:1.7}</style>
    <h1>V4 血栓真实碎裂过程</h1>
    <p>固定相机、Z 轴竖直。读取已接受的 0–3.6 ms 真实计算状态；红色强弱相、蓝灰固定锚点、实际存活键。
    视频时间为逐保存态慢放，右上角标注实际物理时间；没有位置插值或位移放大。</p>
    <h2>原血块区域特写</h2><p>固定毫米坐标窗口。移出该窗口的节点数在左上角显示，完整运动见下一视频。</p>
    <video controls preload="metadata" poster="figures/CLOSEUP_000.png"
    src="videos/VIDEO_V4_FRACTURE_CLOSEUP_PPT_1600.mp4"></video>
    <p><a href="videos/VIDEO_V4_FRACTURE_CLOSEUP_2560.mp4">高清特写 MP4（2560×2048）</a></p>
    <h2>全部碎块完整视野</h2><video controls preload="metadata" src="videos/VIDEO_V4_FRACTURE_FULL_PPT_1600.mp4"></video>
    <p><a href="videos/VIDEO_V4_FRACTURE_FULL_2560.mp4">高清完整视野 MP4（2560×2048）</a></p>
    <h2>真实状态检查图</h2><img src="visual_review/CLOSEUP_CONTACT_SHEET.jpg">
    <img src="visual_review/FULL_CONTACT_SHEET.jpg">
    <p>当前模型有 1927 个材料点、1561 个固定锚点，412 个内部点有初始扰动；边界和锚点仍保留规则排列。
    这是计算模型的碎裂，不代表实验验证。当前机械剪切代理未接入真实声流/CFD 输运。</p>"""
    (OUT/"OPEN_RESULTS.html").write_text(html)
    print("COMPLETE", json.dumps(dict(states=len(states), frames=frames,
                                    wall_s=time.monotonic()-start)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--preview", action="store_true")
    main(parser.parse_args().preview)

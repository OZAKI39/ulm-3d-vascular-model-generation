"""A Z-axis camera orbit of the frozen native state in the loading figure."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time
import io
from itertools import product

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree
import imageio_ffmpeg
from matplotlib import mathtext, rcParams
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.clot_loading_publication_v1.render_loading import C, CASE, V, FONT, sha, write_json, ARROW_MM_PER_UN
from analysis.clot_discretization_three_views_v1.render_three_views import silhouette_axes
from pd_clot.lammps_fragmentation.fragment_analysis import read_dump

SOURCE = ROOT/'results/clot_loading_publication_v1'
OUT = SOURCE/'z_rotation'
WIDTH = HEIGHT = 1920
VIEWPORT = (110, 180, 1700, 1360)  # x, y from top, width, height in pixels
SCALE = 1.05
ELEVATION_DEG = -20.
FPS = 30
FRAMES = 360


def camera(angle, initial):
    theta = np.deg2rad(angle)
    R = np.array([[np.cos(theta), -np.sin(theta), 0.],
                  [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
    outward = R @ initial
    right = np.cross(-outward, [0., 0., 1.]); right /= np.linalg.norm(right)
    up = np.cross(right, -outward)
    return outward, right, up


def project(point, focal, right, up):
    xx, yy, w, h = VIEWPORT
    delta = np.asarray(point)-focal
    return np.array([xx+w/2+(delta@right)*h/(2*SCALE), yy+h/2-(delta@up)*h/(2*SCALE)])


def line_mesh(points, segments):
    return pv.PolyData(points, lines=np.c_[np.full(len(segments), 2), segments].ravel())


def load():
    provenance = json.loads((SOURCE/'PROVENANCE.json').read_text())
    g = dict(np.load(V/'V3_MATERIAL_POINTS_AND_BONDS.npz'))
    tr = np.load(CASE/'accepted_trajectory.npz')
    frame = int(np.flatnonzero(tr['step'] == 540)[0])
    X_m = tr['position_m'][frame]
    table = np.genfromtxt(SOURCE/'selected_frame_nodal_loading.csv', delimiter=',', names=True)
    np.testing.assert_array_equal(table['atom_id'], np.arange(1, len(X_m)+1))
    np.testing.assert_array_equal(np.c_[table['x_m'], table['y_m'], table['z_m']], X_m)
    forces = np.c_[table['load_fx_N'], table['load_fy_N'], table['load_fz_N']]
    for step, cols, a in read_dump(CASE/'trajectory.lammpstrj'):
        if step != 540:
            continue
        a = a[np.argsort(a[:, cols.index('id')])]
        np.testing.assert_array_equal(forces, a[:, [cols.index('v_load_f'+k) for k in 'xyz']])
        break
    else:
        raise AssertionError('Native step 540 missing')
    ev = np.load(CASE/'broken_bonds.npz'); alive = np.ones(len(g['bonds']), bool)
    alive[ev['initial_bond_indices'][ev['events'][:, 1] <= 540]] = False
    assert int(alive.sum()) == provenance['metrics']['native_frame_surviving_bonds']
    arrows = np.array(provenance['rendering']['arrow_atom_ids'])-1
    np.testing.assert_array_equal(np.sort(arrows), np.flatnonzero(table['arrow_shown']))
    assert len(arrows) == 23
    return g, X_m*1e3, forces, g['bonds'][alive], g['bond_types'][alive], arrows, provenance


def build_scene(g, X, forces, bonds, types, arrow_ids, provenance, complete_grid=False):
    bounds = np.array(provenance['rendering']['local_bounds_mm'])
    focal = np.array(provenance['rendering']['camera_focal_mm'])
    original = np.array(provenance['rendering']['camera_outward'])
    azimuth = np.arctan2(original[1], original[0])
    elevation = np.deg2rad(ELEVATION_DEG)
    initial = np.array([np.cos(elevation)*np.cos(azimuth), np.cos(elevation)*np.sin(azimuth), np.sin(elevation)])
    # Preserve the original figure's highlighted links even though the user
    # has requested a lower, more side-on camera angle.
    _, right, up = camera(0, original)
    p = pv.Plotter(off_screen=True, window_size=list(VIEWPORT[2:]), lighting='three lights')
    p.set_background('white'); p.render_window.SetMultiSamples(0); p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=10, occlusion_ratio=.01)
    # Grid planes stay fixed in world coordinates throughout the orbit.
    grid = []
    for plane in range(3):
        other = [k for k in range(3) if k != plane]
        planes = bounds[plane] if complete_grid else [bounds[plane, 0 if original[plane] > 0 else 1]]
        for fixed in planes:
            for varying in other:
                constant = next(k for k in other if k != varying)
                for value in np.linspace(*bounds[constant], 4 if constant == 2 else 5):
                    edge = np.zeros((2, 3)); edge[:, plane] = fixed; edge[:, constant] = value
                    edge[:, varying] = bounds[varying]; grid.append(edge)
    if complete_grid:
        grid = list({tuple(edge.ravel()):edge for edge in grid}.values())
    mesh = line_mesh(np.array(grid).reshape(-1, 3), np.arange(2*len(grid)).reshape(-1, 2))
    p.add_mesh(mesh, name='coordinate_grid', color=C['grid'], lighting=False, line_width=1.1, opacity=.90)
    for axis in range(3):
        other = [k for k in range(3) if k != axis]
        for vals in product(*[bounds[k] for k in other]):
            edge = np.zeros((2, 3)); edge[:, axis] = bounds[axis]
            for k, val in zip(other, vals):
                edge[:, k] = val
            p.add_mesh(pv.Line(*edge), color=C['grid'], lighting=False, line_width=.8, opacity=.55)
    for typ, opacity in [(1, .021), (2, .019), (3, .010)]:
        p.add_mesh(line_mesh(X, bonds[types == typ]), color='#9cabb5', opacity=opacity,
                   line_width=.85, lighting=False)
    # Freeze the very same front-link emphasis used by the reference figure.
    # It remains a single static scene: no arrows or links appear/disappear due
    # to a camera-dependent resampling rule.
    projected = np.c_[(X-focal)@right, (X-focal)@up]; depth = (X-focal)@original
    neighbors = cKDTree(projected).query_ball_point(projected, .055)
    front = np.array([depth[js].max()-depth[i] for i, js in enumerate(neighbors)]) < .045
    length = np.linalg.norm(X[bonds[:, 1]]-X[bonds[:, 0]], axis=1)
    emph = front[bonds].all(axis=1) & (types <= 2) & (length <= .085)
    assert int(emph.sum()) == provenance['rendering']['front_emphasized_real_bonds']
    p.add_mesh(line_mesh(X, bonds[emph]), color=C['bond'], opacity=.42, line_width=1.25, lighting=False)
    n = int(g['clot_count'])
    for color, ids, size, alpha in [(C['anchor'], np.arange(n, len(X)), 3.5, .57),
                                    (C['clot'], np.arange(n), 4.1, .94)]:
        p.add_points(pv.PolyData(X[ids]), color=color, point_size=size, opacity=alpha,
                     render_points_as_spheres=True, ambient=.55, diffuse=.45,
                     specular=.10, specular_power=20, show_scalar_bar=False)
    arrows = []
    for i in arrow_ids:
        length = np.linalg.norm(forces[i])*1e6*ARROW_MM_PER_UN
        arrows.append(pv.Arrow(start=X[i], direction=forces[i], scale=length,
                               tip_length=min(.060/length, .36), tip_radius=.019/length,
                               shaft_radius=.0045/length, shaft_resolution=12, tip_resolution=20))
    p.add_mesh(pv.MultiBlock(arrows).combine(), color=C['force'], smooth_shading=True, ambient=.55, diffuse=.45)
    outward, right, up = camera(0., initial)
    p.camera_position = [(focal+20*outward).tolist(), focal.tolist(), up.tolist()]
    p.enable_parallel_projection(); p.camera.parallel_scale = SCALE
    p.show(auto_close=False)
    return p, bounds, focal, initial


FONTS = {n:ImageFont.truetype(str(FONT), n) for n in [25, 28, 30, 32, 34, 48]}
rcParams['mathtext.fontset'] = 'stix'
FORMULA_BUFFER = io.BytesIO()
mathtext.math_to_image(r'$\mathbf{F}_i^{\mathrm{load}}=A_i\,\tau(\mathbf{x}_i(t),t)\,\mathbf{t}_i^{0}$',
                      FORMULA_BUFFER, prop=FontProperties(size=22), dpi=150, format='png', color=C['text'])
FORMULA_BUFFER.seek(0)
FORMULA = Image.open(FORMULA_BUFFER).convert('RGBA')


def decorate(scene, angle, bounds, focal, right, up):
    im = Image.new('RGB', (WIDTH, HEIGHT), 'white'); im.paste(Image.fromarray(scene).convert('RGB'), VIEWPORT[:2])
    draw = ImageDraw.Draw(im); boxes = []
    def text(xy, value, size=30, color=None, anchor='mm', record=True):
        draw.text(tuple(xy), value, font=FONTS[size], fill=color or C['text'], anchor=anchor)
        if record:
            boxes.append((value, draw.textbbox(tuple(xy), value, font=FONTS[size], anchor=anchor)))
    text((75, 63), 'Localized loading', 48, anchor='lm')
    text((75, 118), 't = 0.270 ms · fixed physical state', 30, color=C['axis'], anchor='lm')
    text((1845, 68), f'Z rotation  {angle:05.1f}°', 30, color=C['axis'], anchor='rm')
    axes = silhouette_axes(bounds, right, up)
    ticks = [[-6.8, -6., -5.2], [-.5, 0., .5], [.6, 1., 1.4]]
    for item in axes:
        k = item['axis']; edge = np.array(item['edge']); base = edge.mean(axis=0)
        normal = np.array(item['outward_screen'])*[1., -1.]
        screen = np.array([project(q, focal, right, up) for q in edge])
        draw.line([tuple(s) for s in screen], fill=C['axis'], width=2)
        for value in ticks[k]:
            point = base.copy(); point[k] = value
            loc = project(point, focal, right, up)
            draw.line([tuple(loc), tuple(loc+normal*6)], fill=C['axis'], width=2)
            text(loc+normal*(85 if k == 2 else 28), f'{value:.1f}' if k == 0 else f'{value:g}', 30)
        if k == 2:
            top = screen[np.argmin(screen[:, 1])]
            text(top+np.array([0, -38]), 'Z (mm)', 34)
        else:
            text(project(base, focal, right, up)+normal*140, 'XY'[k]+' (mm)', 34)
    center = np.array(json.loads((CASE/'request.json').read_text())['load_centers_m'][0])*1e3
    q = project(center, focal, right, up)
    draw.line([(q[0]-10, q[1]), (q[0]+10, q[1])], fill=C['text'], width=3)
    draw.line([(q[0], q[1]-10), (q[0], q[1]+10)], fill=C['text'], width=3)
    # Same proportional force scale as the source; true projected geometry,
    # with a screen-plane 2 microN comparison arrow.
    key_x, key_y = 185, 1510; length_px = .24*VIEWPORT[3]/(2*SCALE)
    end = key_x+length_px
    draw.line([(key_x, key_y), (end-15, key_y)], fill=C['force'], width=5)
    draw.polygon([(end, key_y), (end-18, key_y-9), (end-18, key_y+9)], fill=C['force'])
    text((key_x+length_px/2, key_y+37), '2 µN', 30, C['force'])
    fxy = ((WIDTH-FORMULA.width)//2, 1590-FORMULA.height//2)
    im.paste(FORMULA, fxy, FORMULA)
    boxes.append(('Loading formula', (fxy[0], fxy[1], fxy[0]+FORMULA.width, fxy[1]+FORMULA.height)))
    text((WIDTH/2, 1660), '200 Pa reference budget · local peak 750.26 Pa', 32)
    text((WIDTH/2, 1700), 'Pulse duration 1.1605 ms · Gaussian σ = 0.28 mm', 30, C['axis'])
    items = [('dot', 'Clot nodes', C['clot']), ('dot', 'Fixed anchors', C['anchor']),
             ('line', 'Surviving bonds', C['bond']), ('arrow', 'Load force', C['force']),
             ('plus', 'Fixed load center', C['text'])]
    widths = [43+draw.textlength(label, font=FONTS[28]) for _, label, _ in items]
    x = (WIDTH-sum(widths)-55*(len(items)-1))/2; y = 1780
    for (kind, label, color), width in zip(items, widths):
        if kind == 'dot':
            draw.ellipse((x+6, y-7, x+20, y+7), fill=color)
        elif kind in ['line', 'arrow']:
            draw.line([(x, y), (x+28, y)], fill=color, width=3)
            if kind == 'arrow':
                draw.polygon([(x+30, y), (x+18, y-7), (x+18, y+7)], fill=color)
        else:
            draw.line([(x+3, y), (x+23, y)], fill=color, width=3)
            draw.line([(x+13, y-10), (x+13, y+10)], fill=color, width=3)
        text((x+43, y), label, 28, anchor='lm'); x += width+55
    text((WIDTH/2, 1860), 'View rotation about Z · native nodes, surviving bonds and load vectors', 28, C['axis'])
    clipped = [label for label, b in boxes if min(b[:2]) < 0 or b[2] > WIDTH or b[3] > HEIGHT]
    assert not clipped, (angle, clipped)
    overlaps = []
    for i, (label, a) in enumerate(boxes):
        for other, b in boxes[i+1:]:
            if min(a[2], b[2])-max(a[0], b[0]) > 2 and min(a[3], b[3])-max(a[1], b[1]) > 2:
                overlaps.append([label, other])
    assert not overlaps, (angle, overlaps)
    return im


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--preview', action='store_true'); args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sources = [SOURCE/'Fig_clot_loading.png', SOURCE/'PROVENANCE.json', SOURCE/'selected_frame_nodal_loading.csv',
               V/'V3_MATERIAL_POINTS_AND_BONDS.npz', CASE/'accepted_trajectory.npz', CASE/'broken_bonds.npz',
               CASE/'trajectory.lammpstrj', CASE/'request.json']
    before = {str(path):sha(path) for path in sources}
    g, X, forces, bonds, types, arrows, prov = load()
    p, bounds, focal, initial = build_scene(g, X, forces, bonds, types, arrows, prov)
    corners = np.array(list(product(*bounds)))
    # Inspect every camera's framing and labels before the expensive render.
    blank = np.full((VIEWPORT[3], VIEWPORT[2], 3), 255, np.uint8)
    for angle in range(FRAMES):
        outward, right, up = camera(angle, initial)
        xy = np.array([project(q, focal, right, up) for q in corners])
        assert (xy[:, 0] > VIEWPORT[0]).all() and (xy[:, 0] < VIEWPORT[0]+VIEWPORT[2]).all()
        assert (xy[:, 1] > VIEWPORT[1]).all() and (xy[:, 1] < VIEWPORT[1]+VIEWPORT[3]).all()
        decorate(blank, angle, bounds, focal, right, up)
    print('ALL_360_CAMERA_LAYOUTS_PASS', flush=True)
    selected = np.arange(0, 360, 45) if args.preview else np.arange(FRAMES)
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe(); writer = None
    if not args.preview:
        command = [ffmpeg, '-y', '-hide_banner', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24',
                   '-s:v', f'{WIDTH}x{HEIGHT}', '-r', str(FPS), '-i', '-', '-an', '-c:v', 'libx264',
                   '-preset', 'medium', '-crf', '16', '-pix_fmt', 'yuv420p', '-threads', '4',
                   '-movflags', '+faststart', str(OUT/'Fig_clot_loading_Z_rotation.rendering.mp4')]
        writer = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stills = []; timing = time.perf_counter(); camera_records = []
    try:
        for count, angle in enumerate(selected):
            outward, right, up = camera(float(angle), initial)
            p.camera_position = [(focal+20*outward).tolist(), focal.tolist(), up.tolist()]
            p.camera.parallel_scale = SCALE; p.render()
            scene = p.screenshot(return_img=True)
            im = decorate(scene, float(angle), bounds, focal, right, up)
            if writer:
                writer.stdin.write(np.asarray(im, dtype=np.uint8).tobytes())
            if int(angle) % 45 == 0:
                im.save(OUT/f'view_{int(angle):03d}.png')
                small = im.copy(); small.thumbnail((640, 640)); stills.append(small)
            camera_records.append(dict(frame=int(angle), video_time_s=float(angle/FPS),
                                       orbit_angle_deg=float(angle), focal_mm=focal.tolist(),
                                       eye_mm=(focal+20*outward).tolist(), up=up.tolist(), physical_step=540))
            if (count+1) % (1 if args.preview else 30) == 0:
                print('RENDER', count+1, '/', len(selected), 'elapsed_s', round(time.perf_counter()-timing, 1), flush=True)
        if writer:
            writer.stdin.close(); error = writer.stderr.read().decode(); code = writer.wait()
            assert code == 0, error
            (OUT/'Fig_clot_loading_Z_rotation.rendering.mp4').replace(OUT/'Fig_clot_loading_Z_rotation.mp4')
    finally:
        p.close()
        if writer and writer.poll() is None:
            writer.terminate(); writer.wait()
    sheet = Image.new('RGB', (2560, 1280), 'white')
    for j, im in enumerate(stills):
        sheet.paste(im, ((j % 4)*640, (j//4)*640))
    sheet.save(OUT/'rotation_contact_sheet.jpg', quality=94)
    assert before == {path:sha(path) for path in before}
    metadata = dict(status='PASS', physical_step=540, physical_time_ms=.270, duration_s=12,
                    fps=FPS, frames=FRAMES, resolution=[WIDTH, HEIGHT], rotation_degrees=360,
                    orbit='Camera moves at constant elevation and radius about the line parallel to global Z through the clot viewing center',
                    camera_elevation_deg=ELEVATION_DEG,
                    camera_change='User requested less overhead view; elevation magnitude reduced from55.43 to20deg and magnification increased by1.30/1.05',
                    rotation_axis_point_mm=focal.tolist(), rotation_axis_direction=[0, 0, 1],
                    fixed_world_geometry=True, fixed_load_vectors=True, fixed_arrow_sample=True,
                    arrow_atom_ids=(arrows+1).tolist(), arrow_count=len(arrows),
                    clot_nodes=int(g['clot_count']), anchor_nodes=int(g['anchors'].sum()), surviving_bonds=len(bonds),
                    original_figure_unchanged=True, source_sha256=before, all_360_camera_layouts_pass=True,
                    motion_semantics='Only a viewing orbit, not physical rotation or transient force evolution',
                    loop_semantics='360 frames at one degree increments, last angle359deg; looping advances one degree back to first frame',
                    parallel_scale_mm=SCALE, viewport=list(VIEWPORT), new_solver_runs=0,
                    depth_peeling_maximum=10, depth_peeling_occlusion_ratio=.01,
                    render_elapsed_s=time.perf_counter()-timing, ffmpeg_executable=ffmpeg)
    write_json(OUT/('PREVIEW_PROVENANCE.json' if args.preview else 'PROVENANCE.json'), metadata)
    if not args.preview:
        write_json(OUT/'CAMERA_PATH.json', camera_records)
        (OUT/'README.txt').write_text('''左侧加载示意图绕Z轴环绕视频：1920×1920，30fps，12秒，360°。
保持C1 step540（0.270ms）真实状态，1927个血栓节点、1561个固定锚节点、154367条存活连接，以及原图同一组23个加载力箭头。箭头长度显示比例与原图一致。
相机绕穿过局部视图中心、平行于全局Z方向的轴环绕；按用户要求将倾角从约55.43°减小到20°，改成更接近侧视，并将取景放大约24%。节点和力向量的世界坐标固定，视频不是物理时程。全部节点始终包含在视野中，只可能被正常三维投影遮挡。
共360帧，角度0至359°，接回首帧时仍以1°步进连续循环。字形保持正向，坐标标注沿当前可见外边缘放置；原图的网格平面与节点/键/箭头保持固定，不按视角更换力箭头。
主文件 Fig_clot_loading_Z_rotation.mp4；view_000.png为封面；rotation_contact_sheet.jpg为8个方向的检查图。
源图和所有原始实验文件均未改动，未新增物理求解。PROVENANCE.json和CAMERA_PATH.json记录数据与镜头路径。
''')
        files = [path for path in OUT.iterdir() if path.suffix in ['.mp4', '.png', '.jpg']]
        write_json(OUT/'DELIVERY_SHA256.json', {path.name:dict(bytes=path.stat().st_size, sha256=sha(path)) for path in files})
    print('COMPLETE', OUT, flush=True)


if __name__ == '__main__':
    main()

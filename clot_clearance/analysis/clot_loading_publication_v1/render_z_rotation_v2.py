"""Six-sided grid and clean overlays with smooth, geometrically valid axes."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time
from itertools import product

import numpy as np
from PIL import Image, ImageDraw, ImageColor
import imageio_ffmpeg

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.clot_loading_publication_v1 import render_z_rotation as base

OUT = base.SOURCE/'z_rotation_v2'
FADE_HALF_DEGREES = 12.  # 0.8 second total transition at 30 degrees / second
TICKS = [[-6.8, -6., -5.2], [-.5, 0., .5], [.6, 1., 1.4]]


class AxisTransitions:
    def __init__(self, bounds, focal, initial):
        self.bounds = bounds; self.focal = focal; self.initial = initial
        self.keys = []; self.states = {}
        for angle in range(360):
            _, right, up = base.camera(angle, initial); keys = []
            for item in base.silhouette_axes(bounds, right, up):
                edge = np.array(item['edge']); screen = self.project(edge, right, up)
                delta = screen[1]-screen[0]
                normal = np.array([-delta[1], delta[0]])/np.linalg.norm(delta)
                wanted = np.array(item['outward_screen'])*[1., -1.]
                sign = 1 if normal@wanted > 0 else -1
                key = (item['axis'], *edge.ravel(), sign)
                keys.append(key); self.states[key] = dict(axis=item['axis'], edge=edge, sign=sign)
            self.keys.append(keys)
        self.transitions = [[i for i in range(360) if self.keys[i][axis] != self.keys[i-1][axis]] for axis in range(3)]

    def project(self, points, right, up):
        return np.array([base.project(q, self.focal, right, up) for q in points])

    def groups(self, angle):
        answer = []
        for axis in range(3):
            found = False
            for index in self.transitions[axis]:
                delta = (angle-(index-.5)+180.) % 360.-180.
                if abs(delta) >= FADE_HALF_DEGREES:
                    continue
                # Fade the outgoing edge to zero before changing the edge,
                # then fade the incoming edge in. This avoids two coordinate
                # systems overlapping during the transition. At the exact
                # switch both opacity and its first derivative are zero.
                alpha = .5-.5*np.cos(np.pi*abs(delta)/FADE_HALF_DEGREES)
                key = self.keys[(index-1)%360][axis] if delta < 0 else self.keys[index][axis]
                if alpha > 1e-10:
                    answer.append((key, self.states[key], float(alpha)))
                found = True; break
            if not found:
                key = self.keys[int(np.floor(angle+.5))%360][axis]
                answer.append((key, self.states[key], 1.))
        for axis in range(3):
            assert 0 <= sum(a for _, s, a in answer if s['axis'] == axis) <= 1+1e-12
        return answer

    def labels(self, state, right, up):
        axis = state['axis']; edge = state['edge']; midpoint = edge.mean(axis=0)
        screen = self.project(edge, right, up); d = screen[1]-screen[0]
        normal = state['sign']*np.array([-d[1], d[0]])/np.linalg.norm(d)
        labels = []; ticks = []
        for value in TICKS[axis]:
            point = midpoint.copy(); point[axis] = value
            loc = base.project(point, self.focal, right, up)
            ticks.append((loc, loc+normal*6))
            label = f'{value:.1f}' if axis == 0 else f'{value:g}'
            labels.append((label, loc+normal*(85 if axis == 2 else 28), 30))
        if axis == 2:
            loc = screen[np.argmin(screen[:, 1])]+[0., -38.]
        else:
            loc = base.project(midpoint, self.focal, right, up)+normal*140
        labels.append(('XYZ'[axis]+' (mm)', loc, 34))
        return screen, ticks, labels


def decorate(scene, angle, axes, center):
    im = Image.new('RGBA', (base.WIDTH, base.HEIGHT), 'white')
    im.paste(Image.fromarray(scene).convert('RGBA'), base.VIEWPORT[:2])
    _, right, up = base.camera(angle, axes.initial)
    records = []
    for key, state, alpha in axes.groups(angle):
        layer = Image.new('RGBA', im.size); draw = ImageDraw.Draw(layer)
        line_color = (*ImageColor.getrgb(base.C['axis']), round(255*alpha))
        text_color = (*ImageColor.getrgb(base.C['text']), round(255*alpha))
        edge, ticks, labels = axes.labels(state, right, up)
        draw.line([tuple(q) for q in edge], fill=line_color, width=2)
        for p, q in ticks:
            draw.line([tuple(p), tuple(q)], fill=line_color, width=2)
        for text, position, size in labels:
            draw.text(tuple(position), text, font=base.FONTS[size], fill=text_color, anchor='mm')
            records.append(dict(label=text, xy=position.tolist(), alpha=alpha, state=str(key),
                                bbox=list(draw.textbbox(tuple(position), text, font=base.FONTS[size], anchor='mm'))))
        im.alpha_composite(layer)
    draw = ImageDraw.Draw(im)
    def text(position, label, size, color=None, anchor='mm'):
        draw.text(tuple(position), label, font=base.FONTS[size], fill=color or base.C['text'], anchor=anchor)
        records.append(dict(label=label, xy=list(position), alpha=1., state='fixed_overlay',
                            bbox=list(draw.textbbox(tuple(position), label, font=base.FONTS[size], anchor=anchor))))
    q = base.project(center, axes.focal, right, up)
    draw.line([(q[0]-10, q[1]), (q[0]+10, q[1])], fill=base.C['text'], width=3)
    draw.line([(q[0], q[1]-10), (q[0], q[1]+10)], fill=base.C['text'], width=3)
    key_x, key_y = 185, 1510; length = .24*base.VIEWPORT[3]/(2*base.SCALE); end = key_x+length
    draw.line([(key_x, key_y), (end-15, key_y)], fill=base.C['force'], width=5)
    draw.polygon([(end, key_y), (end-18, key_y-9), (end-18, key_y+9)], fill=base.C['force'])
    text((key_x+length/2, key_y+37), '2 µN', 30, base.C['force'])
    items = [('dot', 'Clot nodes', base.C['clot']), ('dot', 'Fixed anchors', base.C['anchor']),
             ('line', 'Surviving bonds', base.C['bond']), ('arrow', 'Load force', base.C['force']),
             ('plus', 'Fixed load center', base.C['text'])]
    widths = [43+draw.textlength(label, font=base.FONTS[28]) for _, label, _ in items]
    x = (base.WIDTH-sum(widths)-55*(len(items)-1))/2; y = 1780
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
    clipped = [r['label'] for r in records if min(r['bbox'][:2]) < 0 or r['bbox'][2] > base.WIDTH or r['bbox'][3] > base.HEIGHT]
    assert not clipped, (angle, clipped)
    return im.convert('RGB'), records


def layout_checks(axes, center):
    blank = np.full((base.VIEWPORT[3], base.VIEWPORT[2], 3), 255, np.uint8)
    traces = []; overlaps = []; corners = np.array(list(product(*axes.bounds)))
    for angle in range(360):
        _, right, up = base.camera(angle, axes.initial)
        projected = axes.project(corners, right, up)
        x, y, w, h = base.VIEWPORT
        assert (projected[:, 0] > x).all() and (projected[:, 0] < x+w).all()
        assert (projected[:, 1] > y).all() and (projected[:, 1] < y+h).all()
        _, records = decorate(blank, angle, axes, center)
        traces.append({(r['state'], r['label']):r for r in records})
        for i, r in enumerate(records):
            for s in records[i+1:]:
                if min(r['alpha'], s['alpha']) < .15:
                    continue
                a = r['bbox']; b = s['bbox']
                if min(a[2], b[2])-max(a[0], b[0]) > 2 and min(a[3], b[3])-max(a[1], b[1]) > 2:
                    overlaps.append(dict(frame=angle, labels=[r['label'], s['label']]))
    assert not overlaps, overlaps[:12]
    max_opacity_change = 0.; max_visible_position_step = 0.
    for i, curr in enumerate(traces):
        following = traces[(i+1)%360]
        for key in curr.keys() | following.keys():
            a = curr.get(key); b = following.get(key)
            max_opacity_change = max(max_opacity_change, abs((a['alpha'] if a else 0)-(b['alpha'] if b else 0)))
            if a and b and min(a['alpha'], b['alpha']) > .02:
                max_visible_position_step = max(max_visible_position_step, float(np.linalg.norm(np.array(a['xy'])-b['xy'])))
    assert max_opacity_change < .14 and max_visible_position_step < 35
    return dict(status='PASS', checked_frames=360, page_clipping=0, visible_text_overlaps=overlaps,
                transition_centers_deg=[[x-.5 for x in t] for t in axes.transitions],
                transition_duration_s=2*FADE_HALF_DEGREES/30,
                maximum_per_frame_opacity_change=max_opacity_change,
                maximum_per_frame_visible_label_motion_px=max_visible_position_step,
                text_inventory=sorted({r['label'] for trace in traces for r in trace.values()}),
                interpolation='Periodic cosine fade-out then fade-in; opacity and its derivative vanish at the edge switch. Positions remain attached to true box edges; no duplicated coordinate text or fictitious intermediate axes')


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--preview', action='store_true'); args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    old = json.loads((base.OUT/'PROVENANCE.json').read_text())
    protected = [Path(p) for p in old['source_sha256']]+[base.OUT/'Fig_clot_loading_Z_rotation.mp4', base.OUT/'PROVENANCE.json']
    before = {str(p):base.sha(p) for p in protected}
    data = base.load(); g, X, forces, bonds, types, arrow_ids, _ = data
    p, bounds, focal, initial = base.build_scene(*data, complete_grid=True)
    mesh = p.actors['coordinate_grid'].mapper.dataset
    grid_lines = mesh.points[mesh.lines.reshape(-1, 3)[:, 1:]]
    faces = []
    for axis in range(3):
        for side, value in enumerate(bounds[axis]):
            count = int(np.isclose(grid_lines[:, :, axis], value).all(axis=1).sum())
            assert count == (10 if axis == 2 else 9), (axis, side, count)
            faces.append(dict(axis='XYZ'[axis], side=side, position_mm=float(value), grid_segments=count))
    axes = AxisTransitions(bounds, focal, initial)
    center = np.array(json.loads((base.CASE/'request.json').read_text())['load_centers_m'][0])*1e3
    layout = layout_checks(axes, center)
    base.write_json(OUT/'LAYOUT_QA.json', layout)
    print('ALL_360_LAYOUTS_PASS', layout['maximum_per_frame_opacity_change'], layout['maximum_per_frame_visible_label_motion_px'], flush=True)
    selected = sorted(set(range(0, 360, 45)) | {66, 72, 78, 156, 162, 168}) if args.preview else range(360)
    writer = None; ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    if not args.preview:
        writer = subprocess.Popen([ffmpeg, '-y', '-hide_banner', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24',
            '-s:v', '1920x1920', '-r', '30', '-i', '-', '-an', '-c:v', 'libx264', '-preset', 'medium',
            '-crf', '16', '-pix_fmt', 'yuv420p', '-threads', '4', '-movflags', '+faststart',
            str(OUT/'Fig_clot_loading_Z_rotation.rendering.mp4')], stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    camera_records = []; timing = time.perf_counter(); stills = []
    try:
        for count, angle in enumerate(selected):
            outward, right, up = base.camera(angle, initial)
            p.camera_position = [(focal+20*outward).tolist(), focal.tolist(), up.tolist()]
            p.camera.parallel_scale = base.SCALE; p.render()
            im, _ = decorate(p.screenshot(return_img=True), angle, axes, center)
            if writer:
                writer.stdin.write(np.asarray(im, dtype=np.uint8).tobytes())
            if angle % 45 == 0 or args.preview:
                im.save(OUT/f'view_{angle:03d}.png')
            if angle % 45 == 0:
                small = im.copy(); small.thumbnail((640, 640)); stills.append(small)
            camera_records.append(dict(frame=angle, video_time_s=angle/30., orbit_angle_deg=float(angle),
                                       focal_mm=focal.tolist(), eye_mm=(focal+20*outward).tolist(), up=up.tolist(), physical_step=540))
            if args.preview or (count+1) % 30 == 0:
                print('RENDER', count+1, '/', len(selected), 'elapsed_s', round(time.perf_counter()-timing, 1), flush=True)
        if writer:
            writer.stdin.close(); error = writer.stderr.read().decode(); status = writer.wait()
            assert status == 0, error
            (OUT/'Fig_clot_loading_Z_rotation.rendering.mp4').replace(OUT/'Fig_clot_loading_Z_rotation.mp4')
    finally:
        p.close()
        if writer and writer.poll() is None:
            writer.terminate(); writer.wait()
    sheet = Image.new('RGB', (2560, 1280), 'white')
    for i, im in enumerate(stills):
        sheet.paste(im, ((i % 4)*640, (i//4)*640))
    sheet.save(OUT/'rotation_contact_sheet.jpg', quality=94)
    assert before == {name:base.sha(Path(name)) for name in before}
    provenance = {**old, 'source_sha256':before, 'revision':'Six faces, clean labels, smooth axis changes',
        'grid_faces':faces, 'grid_face_count':6, 'grid_line_count':len(grid_lines), 'layout':layout,
        'overlay_text_scope':'Only axis names and tick numbers, shared legend and 2 microN force scale',
        'removed_text':['main title', 'snapshot time', 'rotation angle', 'formula', 'loading parameters', 'footer'],
        'camera_geometry_materials_and_force_sample_match_previous_video':True,
        'previous_video_unchanged':True, 'render_elapsed_s':time.perf_counter()-timing}
    base.write_json(OUT/('PREVIEW_PROVENANCE.json' if args.preview else 'PROVENANCE.json'), provenance)
    if not args.preview:
        base.write_json(OUT/'CAMERA_PATH.json', camera_records)
        (OUT/'README.txt').write_text('''旋转视频 v2：六个面完整网格，精简文字，坐标标注平滑过渡。
沿用上一版20°侧视倾角、取景、配色、1920×1920分辨率、12秒/30fps。血栓、锚节点、连接键和同一组23个力箭头均使用原来0.270ms真实状态。
六个面的网格永久存在；旋转不会留下未绘制网格的面。画面仅保留XYZ坐标轴标题与数字、共享图例、2µN箭头力标尺和加载中心十字标记。删除总标题、时间、转角、公式、参数说明及页脚。
坐标轴换边时使用0.8秒余弦淡出再淡入。换边瞬间透明度与透明度变化率均为零，避免位置突然跳动和两组文字重影。每组坐标始终对应真实包围盒边缘；不会将坐标值移动到虚构的中间轴线上。六面网格与图例、力标尺始终显示。
已检查全部360帧的裁切、文字重叠、标签位置变化与透明度连续性。原视频及源数据保留。
Fig_clot_loading_Z_rotation.mp4为主视频；view_000.png为预览；rotation_contact_sheet.jpg为8视角检查图。
''')
        files = [q for q in OUT.iterdir() if q.suffix in ['.mp4', '.png', '.jpg']]
        base.write_json(OUT/'DELIVERY_SHA256.json', {q.name:dict(bytes=q.stat().st_size, sha256=base.sha(q)) for q in files})
    print('COMPLETE', OUT, flush=True)


if __name__ == '__main__':
    main()

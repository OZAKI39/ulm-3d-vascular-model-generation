"""Re-render candidate 0 streamlines only, reusing the existing visual style.

No CFD, integration, interpolation or particle work is performed. The fixed
camera scale fits the complete 360-degree orbit with an explicit pixel margin.
Render to a staging folder first; --publish installs the verified deliverables.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import time

import imageio_ffmpeg
import numpy as np
import pyvista as pv

import style as shared

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'visualization/candidate_0'
WORK = OUT / 'streamlines_adjustment_20260929'
STAGE = WORK / 'rendered'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load(path):
    return json.loads(Path(path).read_text())


def render(duration_s, margin_px):
    started = time.time()
    active = load(ROOT / 'reports/ACTIVE_FLOW.json')
    prep = load(ROOT / 'visualization/data/PREPARATION.json')
    geometry = load(ROOT / 'inputs/geometry.json')
    assert active['status'] == 'PASS'
    assert prep['source_flow_sha256'] == active['flow_sha256']
    case = ROOT / active['case']
    inputs = {str((case / key).relative_to(ROOT)): value
              for key, value in active['postprocess_data_sha256'].items()}
    inputs.update({str(Path('visualization/data') / key): value
                   for key, value in prep['data_sha256'].items()})
    for path, digest in inputs.items():
        assert sha(ROOT / path) == digest, path

    assert duration_s > 0 and 0 < margin_px < 100
    shared.FRAMES = round(duration_s * shared.FPS)
    assert np.isclose(shared.FRAMES / shared.FPS, duration_s)
    shared.KEYS = sorted(set([0, shared.FRAMES // 8, shared.FRAMES // 4,
                             3 * shared.FRAMES // 8, shared.FRAMES // 2,
                             5 * shared.FRAMES // 8, 3 * shared.FRAMES // 4,
                             7 * shared.FRAMES // 8, shared.FRAMES - 1]))
    shared.SPEED_MAX = float(np.ceil(active['measurements']['velocity_max_m_s'] * 10) * 100)
    T = np.asarray(geometry['pose_transforms']['0']['source_mm_to_candidate_mm'])
    origin = np.asarray(geometry['origin_source_mm'])

    def transform(mesh):
        result = mesh.copy()
        result.points = (np.asarray(result.points, dtype=float) * 1000 + origin) @ T[:3, :3].T + T[:3, 3]
        return result

    wall = transform(pv.read(case / 'postprocess/pressure_surface_si.vtp'))
    lines = transform(pv.read(ROOT / 'visualization/data/streamlines_si.vtp'))
    count = len(prep['selected_streamlines'])
    assert count == len(np.unique(lines.cell_data['Line_id']))
    ports = []
    for name, spec in geometry['ports'].items():
        cap = pv.read(ROOT / 'mesh/SV_MESH/mesh-surfaces' / (spec['role'] + '.vtp'))
        center = (np.asarray(cap.center) * 1000 + origin) @ T[:3, :3].T + T[:3, 3]
        ports.append(('Inlet' if name == 'I1' else 'Outlet 0' + name[-1], center))

    center = shared.fitted_z_center(wall.points)
    baseline = shared.Turntable(wall.points, [0, 0, 1], center)
    # Include line vertices as well as the transparent wall. These analytical
    # bounds hold at every azimuth, not just the frames selected for inspection.
    q = np.vstack([wall.points, lines.points]) - center
    radial = np.linalg.norm(q[:, :2], axis=1)
    phi = np.deg2rad(shared.ELEVATION_DEG)
    half_width = shared.SIZE[0] * .90 / 2
    half_height = shared.SIZE[1] * shared.VIEW_HEIGHT / 2
    horizontal_scale = radial.max() / shared.ASPECT / (1 - margin_px / half_width)
    vertical_scale = (np.abs(q[:, 2]) * np.cos(phi) + radial * np.sin(phi)).max() / (1 - margin_px / half_height)
    fitted_scale = float(max(horizontal_scale, vertical_scale))
    magnification = baseline.scale / fitted_scale
    orbit = shared.Turntable(wall.points, [0, 0, 1], center, magnification=magnification)
    assert np.isclose(orbit.scale, fitted_scale)
    bounds = np.asarray([row['projected_bounds'] for row in orbit.trace])
    assert np.max(np.abs(bounds[:, :2])) <= 1 - margin_px / half_width + 1e-12
    assert np.max(np.abs(bounds[:, 2:])) <= 1 - margin_px / half_height + 1e-12
    print(json.dumps(dict(stage='CAMERA_FIT',magnification=magnification,
                          parallel_scale_mm=orbit.scale, duration_s=duration_s,
                          frames=shared.FRAMES, margin_px=margin_px)), flush=True)
    for folder in ['figures', 'animations', 'inspection']:
        (STAGE / folder).mkdir(parents=True, exist_ok=True)

    view = shared.View(wall, lines, orbit, size=shared.SIZE,
                       kind='full_vessel', count=count, ports=ports)
    annotations = view.annotations
    # This renderer uses the unmodified label planner and overlay. No hidden
    # ports, text resizing, color changes or viewport changes are necessary.
    gpu = [line for line in view.p.render_window.ReportCapabilities().splitlines()
           if 'OpenGL vendor' in line or 'OpenGL renderer' in line]
    dump(WORK / 'CAMERA_FIT.json', dict(
        old_parallel_scale_mm=baseline.scale, new_parallel_scale_mm=orbit.scale,
        magnification_vs_previous=magnification,
        edge_margin_at_1080p_px=margin_px, fixed_center_mm=center.tolist(),
        continuous_orbit_bound=True, bounds_of_all_rendered_frames=bounds.tolist(),
        source_flow_sha256=active['flow_sha256'], annotations=annotations.record))
    for i in shared.KEYS:
        view.frame(i).save(STAGE / 'inspection' / f'streamlines_{i:04d}.png')
    print('PREVIEWS_READY', flush=True)
    path = STAGE / 'animations/streamlines_rotation.mp4'
    tmp = path.with_suffix('.tmp.mp4')
    encoder = imageio_ffmpeg.write_frames(
        str(tmp), size=shared.SIZE, fps=shared.FPS, codec='libx264', quality=8,
        pix_fmt_in='rgb24', pix_fmt_out='yuv420p', macro_block_size=8,
        output_params=['-preset', 'medium', '-threads', '4', '-movflags', '+faststart'])
    encoder.send(None)
    try:
        for i in range(shared.FRAMES):
            encoder.send(np.ascontiguousarray(view.frame(i)))
            if i % 108 == 0:
                print('Rendered', i, '/', shared.FRAMES, flush=True)
    finally:
        encoder.close()
        view.close()
    reader = imageio_ffmpeg.read_frames(str(tmp))
    meta = next(reader)
    decoded = sum(1 for _ in reader)
    assert meta['size'] == shared.SIZE and meta['fps'] == shared.FPS
    assert decoded == shared.FRAMES and np.isclose(meta['duration'], duration_s, atol=.05)
    tmp.replace(path)

    view = shared.View(wall, lines, orbit, size=(3840, 2160), kind='full_vessel',
                       count=count, ports=ports, annotations=annotations)
    image = view.frame(0, still=True)
    image.save(STAGE / 'figures/streamlines_overview_4k.png', dpi=(300, 300))
    image.save(STAGE / 'figures/streamlines_overview.pdf', resolution=300)
    view.close()
    for relative, digest in inputs.items():
        assert sha(ROOT / relative) == digest, relative
    outputs = {str(p.relative_to(STAGE)): sha(p) for p in sorted(STAGE.rglob('*')) if p.is_file()}
    record = dict(
        status='PASS', source_flow_sha256=active['flow_sha256'],
        duration_s=duration_s, fps=shared.FPS, frames=decoded, size=list(shared.SIZE),
        rotation_deg_s=360 / duration_s, old_duration_s=18,
        magnification_vs_previous=magnification, margin_px=margin_px,
        all_azimuths_fit=True, annotation_validation=annotations.record,
        OpenGL=gpu, inputs_sha256=inputs, outputs_sha256=outputs,
        raw_values_unchanged=True, style_code_sha256=sha(ROOT / 'visualization/style.py'),
        renderer_sha256=sha(Path(__file__)), elapsed_s=time.time() - started,
        user_scope='Largest fixed full-orbit fit without cropping; half previous rotation speed',
        scope='Presentation only. Existing streamlines reused; no CFD or microbubble work.')
    dump(WORK / 'RENDER_VALIDATION.json', record)
    print('RENDER_PASS', json.dumps({k: record[k] for k in ['frames','duration_s','magnification_vs_previous','elapsed_s']}), flush=True)


def publish():
    record = load(WORK / 'RENDER_VALIDATION.json')
    assert record['status'] == 'PASS'
    assert load(ROOT / 'reports/ACTIVE_FLOW.json')['flow_sha256'] == record['source_flow_sha256']
    for relative, digest in record['outputs_sha256'].items():
        assert sha(STAGE / relative) == digest, relative
    for relative, digest in record['inputs_sha256'].items():
        assert sha(ROOT / relative) == digest, relative
    # Keep old inspection frames in the existing backup; avoid mixed timebases.
    backup = WORK / 'previous_inspection'
    backup.mkdir(exist_ok=True)
    for path in (OUT / 'inspection').glob('streamlines_*.png'):
        if not (backup / path.name).exists():
            shutil.copy2(path, backup / path.name)
        path.unlink()
    for relative in record['outputs_sha256']:
        dest = OUT / relative
        tmp = dest.with_name(dest.name + '.replacement')
        shutil.copy2(STAGE / relative, tmp)
        tmp.replace(dest)
    manifest = load(OUT / 'MEDIA_VALIDATION.json')
    manifest.pop('frames', None)  # The five movies now have different durations.
    manifest['frame_counts'] = 'See each media record; streamlines updated independently.'
    for row in manifest['media']:
        row.setdefault('fps', manifest['fps'])
        row.setdefault('duration_s', row['frames'] / row['fps'])
        if row['kind'] == 'streamlines':
            row.update(sha256=sha(OUT / row['file']), frames=record['frames'],
                       duration_s=record['duration_s'], rotation_deg_s=record['rotation_deg_s'],
                       magnification_vs_previous=record['magnification_vs_previous'],
                       full_orbit_in_frame=True, OpenGL=record['OpenGL'],
                       validation='streamlines_adjustment_20260929/RENDER_VALIDATION.json')
    dump(OUT / 'MEDIA_VALIDATION.json', manifest)
    html_path = OUT / 'OPEN_RESULTS.html'
    html = html_path.read_text()
    import re
    html = re.sub(r'animations/streamlines_rotation\.mp4(?:\?[^"<>]*)?',
                  'animations/streamlines_rotation.mp4?v=' + sha(OUT / 'animations/streamlines_rotation.mp4')[:12], html)
    html_path.write_text(html)
    print('PUBLISHED', OUT, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--duration-s', type=float, default=36.)
    parser.add_argument('--margin-px', type=float, default=8.)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    if args.publish:
        publish()
    else:
        render(args.duration_s, args.margin_px)

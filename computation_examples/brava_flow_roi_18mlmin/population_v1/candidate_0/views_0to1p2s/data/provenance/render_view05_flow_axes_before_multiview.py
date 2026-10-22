"""Fixed view 05, flow-style axes, recorded candidate-0 particles only.

This renderer never imports a solver. Speed filtering changes display visibility
only; original accepted trajectory endpoints and their velocities remain intact.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import os

os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW', 'vtkEGLRenderWindow')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
OUT = ROOT/'candidate_0/view05_flow_axes'
ANALYSIS = ROOT/'analysis/realization_0000_contactfix'
SIZE = (1920, 1080)
FPS = 24


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def project(plot, points):
    matrix = plot.camera.GetCompositeProjectionTransformMatrix(plot.renderer.GetTiledAspectRatio(), -1, 1)
    matrix = np.array([[matrix.GetElement(i, j) for j in range(4)] for i in range(4)])
    h = np.column_stack((points, np.ones(len(points)))) @ matrix.T
    ndc = h[:, :2]/h[:, 3, None]
    return np.column_stack(((ndc[:, 0]+1)*960, (1-ndc[:, 1])*540))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview', action='store_true')
    parser.add_argument('--end-time', type=float, default=6.)
    args = parser.parse_args()
    for name in ['figures', 'animations', 'data']:
        (OUT/name).mkdir(parents=True, exist_ok=True)
    identity_path = ROOT/'runs/realization_0000_contactfix/IDENTITY.json'
    identity = json.loads(identity_path.read_text())
    catalog_path = ANALYSIS/'trajectory_catalog.csv'
    rows = list(csv.DictReader(catalog_path.open()))
    risk_path = ROOT/'cfd_wall_audit/acceleration_overlap.csv'
    risk_ids = {int(r['particle_id']) for r in csv.DictReader(risk_path.open())
                if r['window_any_eta_gt_0p5'] == 'True'}
    failure_ids = {int(r['particle_id']) for r in rows if r['status'] == 'SOLVER_FAILURE'}
    ids = np.array([int(r['particle_id']) for r in rows])
    prior = np.isin(ids, sorted(risk_ids | failure_ids))
    audit = json.loads((ROOT/'cfd_wall_audit/AUDIT_IDENTITY.json').read_text())
    assert audit['sha256']['flow'] == identity['flow_sha256']
    snapshot_root = ROOT/'candidate_0/static_t6s_white_frame'
    old = json.loads((snapshot_root/'data/manifest.json').read_text())
    reference = old['views'][4]
    assert reference['name'] == 'view_05_az240_el25.0_t6p0s'
    direction = np.array(reference['camera_position'])-reference['camera_focal_point']
    direction /= np.linalg.norm(direction)
    config = json.loads((ROOT/'data/environment.json').read_text())
    if not Path(config['flow']).exists():
        config = json.loads((ROOT/'data/environment.local.json').read_text())
    assert digest(config['flow']) == identity['flow_sha256']
    protected_paths = [identity_path, catalog_path, risk_path, Path(config['flow']),
                       ROOT/'candidate_0/animations/microbubble_physical_time.mp4',
                       ROOT/'candidate_0/data/render_manifest.json',
                       snapshot_root/'data/manifest.json',
                       snapshot_root/'figures'/f"{reference['name']}.png",
                       ROOT/'candidate_15/animations/microbubble_physical_time.mp4',
                       ROOT/'candidate_15/data/render_manifest.json']
    protected = {str(p): digest(p) for p in protected_paths}
    geometry = json.loads((BASE/'inputs/geometry.json').read_text())
    transform = np.array(geometry['pose_transforms']['0']['source_mm_to_candidate_mm'])
    origin = np.array(geometry['origin_source_mm'])

    def pose(x):
        return (np.asarray(x)*1000+origin) @ transform[:3, :3].T + transform[:3, 3]

    surface = pv.read(config['flow']).extract_surface(algorithm='dataset_surface')
    surface.points = pose(surface.points)
    center = np.array(surface.center)
    radius = np.linalg.norm(np.ptp(surface.points, axis=0))/2
    ports = []
    port_records = {}
    for role in ['INLET', 'OUTLET_01', 'OUTLET_02', 'OUTLET_03']:
        cap_path = BASE/'mesh/SV_MESH/mesh-surfaces'/f'{role}.vtp'
        cap = pv.read(cap_path)
        name = 'Inlet' if role == 'INLET' else 'Outlet '+str(int(role[-2:]))
        ports.append((name, pose(np.array(cap.center))))
        port_records[role] = dict(center_mm=ports[-1][1].tolist(), sha256=digest(cap_path))
        if role == 'INLET':
            inlet_area = cap.area
    q = identity['population_contract']['source']['Q_in_m3_s']
    dt_ms = identity['population_contract']['nominal_dt_s']*1000

    # Reuse the current flow visualization's actual CoordinateGrid and overlay.
    # Only its camera projector/available pixel rectangle differ for this fixed
    # view; line colors, widths, opacity, physical ticks and font sizes are kept.
    style_path = BASE/'visualization/style.py'
    spec = importlib.util.spec_from_file_location('flow_coordinate_style', style_path)
    flow_style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(flow_style)
    flow_style.FRAMES = 1
    plot = pv.Plotter(off_screen=True, window_size=SIZE)
    plot.set_background('#000000')
    plot.add_mesh(surface, color='#b2becf', opacity=.16, smooth_shading=True)
    plot.camera.position = center+direction*radius*4
    plot.camera.focal_point = center
    plot.camera.up = (0, 0, 1)
    plot.camera.parallel_projection = True
    plot.camera.parallel_scale = radius
    plot.show(auto_close=False)
    axes = flow_style.CoordinateGrid(plot, surface.bounds, center, 'full_vessel', SIZE)
    forward = -direction
    right = np.cross(forward, [0, 0, 1]); right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    # Fit the vessel, not a bounding cage. Grid planes may continue outside the
    # viewport as in the flow renderer; all vessel points stay inside safe area.
    safe = np.array([[125., 265.], [1795., 975.]])
    for _ in range(3):
        xy = project(plot, surface.points)
        plot.camera.parallel_scale *= np.max(np.ptp(xy, axis=0)/(safe[1]-safe[0]))
        xy = project(plot, surface.points)
        delta = safe.mean(axis=0)-(xy.min(0)+xy.max(0))/2
        shift = (-delta[0]*right+delta[1]*up)*(2*plot.camera.parallel_scale/1080)
        plot.camera.position = np.array(plot.camera.position)+shift
        plot.camera.focal_point = np.array(plot.camera.focal_point)+shift
    xy = project(plot, surface.points)
    vessel_bounds = np.array([xy.min(0), xy.max(0)])
    assert np.all(vessel_bounds[0] >= safe[0]-1e-7)
    assert np.all(vessel_bounds[1] <= safe[1]+1e-7)

    class FixedAnnotations(flow_style.AnnotationPlan):
        def project(self, points, i):
            return project(plot, np.asarray(points))

        @staticmethod
        def in_view(centers, width, height):
            return ((centers[:, :, 0]-width/2 > 26) & (centers[:, :, 0]+width/2 < 1894)
                    & (centers[:, :, 1]-height/2 > 204) & (centers[:, :, 1]+height/2 < 1054))

    annotations = FixedAnnotations(None, surface, axes, ports=None)
    annotations.export(OUT/'data/coordinate_annotations.json')
    # Preserve population port typography; place neutral arrows at real caps.
    titlefont = ImageFont.truetype(str(ROOT/'assets/NotoSansCJKsc-Regular.otf'), 30)
    portfont = ImageFont.truetype(str(ROOT/'assets/NotoSansCJKsc-Regular.otf'), 23)
    occupied = [r['bbox_px'] for r in json.loads((OUT/'data/coordinate_annotations.json').read_text())['frames'][0]['labels'] if r['opacity'] > 0]
    tree = cKDTree(xy)
    draw = ImageDraw.Draw(Image.new('RGB', SIZE))
    port_labels = []
    for name, position in ports:
        cap = project(plot, np.array([position]))[0]
        candidates = []
        for distance in [85, 110, 140, 180, 230]:
            for angle in np.arange(0, 360, 15):
                p = cap+distance*np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))])
                bb = draw.textbbox(tuple(p), name, font=portfont, anchor='mm')
                if bb[0] < 25 or bb[2] > 1895 or bb[1] < 207 or bb[3] > 1055:
                    continue
                if any(bb[0]-8 < b[2] and bb[2]+8 > b[0] and bb[1]-8 < b[3] and bb[3]+8 > b[1] for b in occupied):
                    continue
                sample = np.array(np.meshgrid(np.linspace(bb[0]-8, bb[2]+8, 20), np.linspace(bb[1]-8, bb[3]+8, 6))).reshape(2, -1).T
                clearance = tree.query(sample)[0].min()
                if clearance < 12:
                    continue
                tail = np.clip(cap, np.array(bb[:2])-6, np.array(bb[2:])+6)
                candidates.append((distance, p, bb, tail))
        assert candidates, name
        _, p, bb, tail = min(candidates, key=lambda c: c[0])
        occupied.append(bb)
        port_labels.append(dict(text=name, center_px=p.tolist(), bbox_px=list(bb), cap_px=cap.tolist(), tail_px=tail.tolist()))

    def overlay(pixels, t):
        img = annotations.overlay(Image.fromarray(pixels), 0, still=True)
        draw = ImageDraw.Draw(img)
        for row in port_labels:
            cap = np.array(row['cap_px']); tail = np.array(row['tail_px'])
            draw.text(tuple(row['center_px']), row['text'], font=portfont, anchor='mm', fill='white')
            draw.line([tuple(tail), tuple(cap)], fill='#dce2ed', width=2)
            v = tail-cap; v /= np.linalg.norm(v); side = np.array([-v[1], v[0]])
            draw.polygon([tuple(cap), tuple(cap+12*v+4*side), tuple(cap+12*v-4*side)], fill='#dce2ed')
        for i, line in enumerate(['Continuous microbubble perfusion', f'Physical time: {t:.3f} s',
                                  f'dt = {dt_ms:g} ms', f'Inlet mean speed: {q/inlet_area:.3f} m/s ({q*60e6:g} mL/min)']):
            draw.text((40, 22+40*i), line, font=titlefont, fill='white')
        return img

    # Snapshot first validates camera/layout independently of animation sampling.
    snap = np.load(snapshot_root/'data/snapshot_particles.npz')
    cloud = pv.PolyData(pose(snap['position_m'][snap['displayed']]))
    actor = plot.add_mesh(cloud, color='#3399ff', point_size=9.1, render_points_as_spheres=True, lighting=True)
    plot.render()
    overlay(plot.screenshot(return_img=True), 6.).save(OUT/'figures/view05_t6p0s.png', dpi=(300, 300))
    gpu = [r for r in plot.ren_win.ReportCapabilities().splitlines() if any(s in r for s in ['vendor string', 'renderer string', 'version string'])]
    assert any('NVIDIA' in r for r in gpu), gpu
    camera_error = float(np.linalg.norm((np.array(plot.camera.position)-plot.camera.focal_point)/np.linalg.norm(np.array(plot.camera.position)-plot.camera.focal_point)-direction))
    assert camera_error < 1e-12
    layout = dict(reference_view=reference['name'], fixed_camera=True, direction_error=camera_error,
                  camera_position=list(plot.camera.position), camera_focal_point=list(plot.camera.focal_point),
                  parallel_scale=plot.camera.parallel_scale, zoom_relative_to_reference=reference['parallel_scale']/plot.camera.parallel_scale,
                  vessel_bounds_px=vessel_bounds.tolist(), safe_vessel_bounds_px=safe.tolist(),
                  ports=port_records, port_labels=port_labels, white_frame=False,
                  grid_source=str(style_path), grid_source_sha256=digest(style_path),
                  grid_origin_mm=axes.origin.tolist(), tick_values_mm=[x.tolist() for x in axes.ticks],
                  render_backend=gpu)
    dump(OUT/'data/layout.json', layout)
    print('LAYOUT', json.dumps(dict(zoom=layout['zoom_relative_to_reference'], vessel_bounds_px=layout['vessel_bounds_px'])), flush=True)
    if args.preview:
        plot.close()
        return

    old_animation = json.loads((ROOT/'candidate_0/data/render_manifest.json').read_text())
    seconds_per_frame = old_animation['physical_time_end_s']/(old_animation['frames']-1)
    assert 0 < args.end_time <= old_animation['physical_time_end_s']
    times = np.linspace(0, args.end_time, round(args.end_time/seconds_per_frame)+1)
    nframes = len(times)
    positions = np.full((nframes, len(rows), 3), np.nan)
    speeds = np.full((nframes, len(rows)), np.nan)
    durations = np.zeros_like(speeds)
    previous_hashes = json.loads((ROOT/'candidate_0/data/animation_track_hashes.json').read_text())
    hashes = {}
    maximum_speed_error = 0.
    for i, row in enumerate(rows):
        path = ANALYSIS/'tracks'/f'mb_{ids[i]:06d}.npz'
        hashes[path.name] = digest(path)
        assert hashes[path.name] == previous_hashes[path.name]
        with np.load(path) as track:
            t, x, v = track['physical_time_s'], track['x_m'], track['velocity_m_s']
            assert np.all(np.diff(t) > 0)
            live = (times >= t[0]) & (times < t[-1])
            live_t = times[live]
            for component in range(3):
                positions[live, i, component] = np.interp(live_t, t, x[:, component])
            # Stored v[k+1] is the accepted held velocity on [t[k],t[k+1]].
            speed = np.linalg.norm(v[1:], axis=1)
            independent = np.linalg.norm(np.diff(x, axis=0), axis=1)/np.diff(t)
            maximum_speed_error = max(maximum_speed_error, float(np.max(abs(speed-independent))))
            integral = np.r_[0., np.cumsum(speed*np.diff(t))]
            lo = np.maximum(t[0], live_t-.05); hi = np.minimum(t[-1], live_t+.05)
            durations[live, i] = hi-lo
            speeds[live, i] = (np.interp(hi, t, integral)-np.interp(lo, t, integral))/(hi-lo)
        if i % 250 == 0:
            print('READ_TRACK', i, len(rows), flush=True)
    assert maximum_speed_error < 1e-8, maximum_speed_error
    alive = np.isfinite(positions[:, :, 0])
    shown = alive & ~prior[None, :] & (speeds >= .002)
    assert not np.any(shown[:, prior])
    t6_match = None
    if abs(times[-1]-6.) < 1e-12:
        t6_match = bool(np.array_equal(ids[shown[-1]], np.array(old['displayed_particle_ids'])))
        assert t6_match
    np.savez_compressed(OUT/'data/animation_samples.npz', frame_time_s=times, particle_id=ids,
                        position_m=positions, mean_speed_m_s=speeds, averaging_duration_s=durations,
                        displayed=shown, prior_exclusion=prior)
    dump(OUT/'data/track_hashes.json', hashes)
    with (OUT/'data/frame_summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['frame', 'video_time_s', 'physical_time_s', 'alive', 'prior_hidden', 'slow_hidden_additional', 'displayed'])
        for j, t in enumerate(times):
            writer.writerow([j, j/FPS, t, int(alive[j].sum()), int((alive[j] & prior).sum()),
                             int((alive[j] & ~prior & ~shown[j]).sum()), int(shown[j].sum())])
    # Keep one point per source ID; NaN points represent hidden/inactive particles.
    plot.remove_actor(actor)
    cloud = pv.PolyData(np.zeros((len(ids), 3)))
    actor = plot.add_mesh(cloud, color='#3399ff', point_size=9.1, render_points_as_spheres=True, lighting=True)
    import imageio.v2 as imageio
    os.environ['IMAGEIO_FFMPEG_EXE'] = '/usr/bin/ffmpeg'
    movie = OUT/'animations/microbubble_view05_physical_time.mp4'
    writer = imageio.get_writer(movie, fps=FPS, codec='libx264', quality=8, macro_block_size=1, ffmpeg_params=['-threads', '2'])
    frames_mm = pose(positions)
    frames_mm[~shown] = np.nan
    for j, t in enumerate(times):
        cloud.points = frames_mm[j]
        plot.render()
        frame = overlay(plot.screenshot(return_img=True), float(t))
        writer.append_data(np.asarray(frame))
        if j in [0, nframes//4, nframes//2, nframes-1]:
            frame.save(OUT/'figures'/f'frame_{j:04d}.png')
        if j % 96 == 0:
            print('FRAME', j, nframes, flush=True)
    writer.close(); plot.close()
    assert all(digest(p) == value for p, value in protected.items())
    dump(OUT/'data/protected_files_sha256.json', protected)
    record = dict(candidate_id=0, run='realization_0000_contactfix', flow_sha256=identity['flow_sha256'],
                  renderer_sha256=digest(__file__), layout=layout, frames=nframes, fps=FPS,
                  physical_time_s=[float(times[0]), float(times[-1])], video_duration_s=nframes/FPS,
                  physical_seconds_per_video_second=(times[-1]-times[0])/((nframes-1)/FPS),
                  preserved_slowdown_relative_to_original=.2, source_particle_count=len(ids),
                  prior_hidden_ids=ids[prior].tolist(), failure_ids=sorted(failure_ids), risk_ids=sorted(risk_ids),
                  low_speed_filter=dict(threshold_mm_s=2., window_offsets_s=[-.05, .05],
                    definition='Time-weighted mean of true 3-D speed magnitude on accepted segments',
                    lifetime_edges='Use available lifetime overlap; no zero filling or extrapolation',
                    scope='Per-frame display filter only; particles may appear/disappear on threshold crossings'),
                  shown_anytime_count=int(shown.any(0).sum()), shown_at_final_frame=int(shown[-1].sum()),
                  t6_particle_ids_match_static_snapshot=t6_match,
                  maximum_velocity_vs_displacement_speed_error_m_s=maximum_speed_error,
                  protected_files_unchanged=True, text_language='English', background='#000000',
                  bubble_color='#3399ff', marker_size_px=9.1, vessel_opacity=.16,
                  movie_sha256=digest(movie), simulation_or_CFD_started=False,
                  limitation='Display filtering does not certify remaining trajectories physically accurate')
    dump(OUT/'data/render_manifest.json', record)
    print('COMPLETE', json.dumps({k: record[k] for k in ['frames', 'video_duration_s', 'shown_anytime_count', 'shown_at_final_frame']}), flush=True)


if __name__ == '__main__':
    main()

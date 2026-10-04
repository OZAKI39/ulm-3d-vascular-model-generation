"""Candidate 0 snapshots from recorded trajectories; no flow/particle solves."""
from pathlib import Path
import csv
import hashlib
import json
import os
import math

os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW', 'vtkEGLRenderWindow')
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree
from spatial_axes_reference import add_spatial_reference_axes, axis_annotations
from style_reference import projected_bounds

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
RUN = 'realization_0000_contactfix'
ANALYSIS = ROOT / 'analysis' / RUN
OUT = ROOT / 'candidate_0/static_t6s_white_frame'
for directory in [OUT, OUT/'figures', OUT/'data']:
    directory.mkdir(parents=True, exist_ok=True)
TIME = 6.0
HALF_WINDOW = .05
MIN_SPEED_M_S = .002

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

catalog = list(csv.DictReader((ANALYSIS/'trajectory_catalog.csv').open()))
risk_path = ROOT/'cfd_wall_audit/acceleration_overlap.csv'
risks = {int(row['particle_id']) for row in csv.DictReader(risk_path.open())
         if row['window_any_eta_gt_0p5'] == 'True'}
failures = {int(row['particle_id']) for row in catalog if row['status'] == 'SOLVER_FAILURE'}
identity = json.loads((ROOT/'runs'/RUN/'IDENTITY.json').read_text())
audit = json.loads((ROOT/'cfd_wall_audit/AUDIT_IDENTITY.json').read_text())
assert audit['sha256']['flow'] == identity['flow_sha256']
source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in [
    ANALYSIS/'trajectory_catalog.csv', risk_path, ROOT/'runs'/RUN/'IDENTITY.json']}
previous_track_hashes = json.loads((ROOT/'candidate_0/data/animation_track_hashes.json').read_text())
rows, active_ids, active_x, active_mean_speed, visible = [], [], [], [], []
track_hashes = {}
maximum_speed_check_error = 0.
maximum_position_check_error = 0.
partial_windows = []
for item in catalog:
    pid = int(item['particle_id'])
    alive = float(item['birth_time_s']) <= TIME < float(item['physical_time_s'])
    reasons = []
    if pid in failures:
        reasons.append('SOLVER_FAILURE')
    if pid in risks:
        reasons.append('AUDITED_SPEEDUP_WINDOW_NUMERICAL_RISK')
    row = dict(particle_id=pid, active_at_6s=alive, displayed=False,
               prior_solver_failure=pid in failures, prior_numerical_risk=pid in risks,
               slow_below_2_mm_s=False, mean_speed_mm_s='', instantaneous_speed_mm_s='',
               averaging_start_s='', averaging_end_s='', averaging_duration_s='',
               full_0p1s_window=False, x_m='', y_m='', z_m='', hidden_reason='')
    if alive:
        path = ANALYSIS/'tracks'/f'mb_{pid:06d}.npz'
        track_hashes[path.name] = digest(path)
        assert track_hashes[path.name] == previous_track_hashes[path.name]
        with np.load(path) as track:
            t, x, v = track['physical_time_s'], track['x_m'], track['velocity_m_s']
            assert np.all(np.diff(t) > 0)
            # v[i+1] is the held velocity for [t[i], t[i+1]]. Average its
            # magnitude by actual overlapping duration, not snapshot count.
            weights = np.maximum(0., np.minimum(t[1:], TIME+HALF_WINDOW)
                                  - np.maximum(t[:-1], TIME-HALF_WINDOW))
            use = weights > 0
            duration = float(weights.sum())
            mean_speed = float(np.dot(weights[use], np.linalg.norm(v[1:][use], axis=1))/duration)
            path_speed = np.linalg.norm(np.diff(x, axis=0)[use], axis=1)/np.diff(t)[use]
            independent = float(np.dot(weights[use], path_speed)/duration)
            maximum_speed_check_error = max(maximum_speed_check_error, abs(mean_speed-independent))
            segment = int(np.searchsorted(t, TIME, side='right')-1)
            position = np.array([np.interp(TIME, t, x[:,k]) for k in range(3)])
            held_position = x[segment] + (TIME-t[segment])*v[segment+1]
            maximum_position_check_error = max(maximum_position_check_error,
                                                float(np.linalg.norm(position-held_position)))
            slow = mean_speed < MIN_SPEED_M_S
            if slow:
                reasons.append('WINDOW_MEAN_SPEED_LT_2_MM_S')
            full_window = abs(duration-.1) < 1e-10
            if not full_window:
                partial_windows.append(pid)
            row.update(displayed=not reasons, slow_below_2_mm_s=slow,
                       mean_speed_mm_s=mean_speed*1000,
                       instantaneous_speed_mm_s=float(np.linalg.norm(v[segment+1]))*1000,
                       averaging_start_s=max(float(t[0]), TIME-HALF_WINDOW),
                       averaging_end_s=min(float(t[-1]), TIME+HALF_WINDOW),
                       averaging_duration_s=duration, full_0p1s_window=full_window,
                       x_m=position[0], y_m=position[1], z_m=position[2])
            active_ids.append(pid); active_x.append(position)
            active_mean_speed.append(mean_speed); visible.append(not reasons)
    else:
        reasons.append('NOT_ACTIVE_AT_SNAPSHOT')
    row['hidden_reason'] = ';'.join(reasons)
    rows.append(row)

assert maximum_speed_check_error < 1e-11
assert maximum_position_check_error < 1e-12
active_ids = np.asarray(active_ids)
active_x = np.asarray(active_x)
visible = np.asarray(visible)
shown_ids = active_ids[visible]
assert not (set(shown_ids) & (risks | failures))
with (OUT/'data/display_selection.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
np.savez_compressed(OUT/'data/snapshot_particles.npz', physical_time_s=TIME,
                    particle_id=active_ids, position_m=active_x,
                    mean_speed_m_s=np.asarray(active_mean_speed), displayed=visible)

geometry = json.loads((BASE/'inputs/geometry.json').read_text())
transform = np.array(geometry['pose_transforms']['0']['source_mm_to_candidate_mm'])
origin = np.array(geometry['origin_source_mm'])
def pose(x):
    return (np.asarray(x)*1000+origin) @ transform[:3,:3].T + transform[:3,3]

config = json.loads((ROOT/'data/environment.json').read_text())
if not Path(config['flow']).exists():
    config = json.loads((ROOT/'data/environment.local.json').read_text())
assert digest(Path(config['flow'])) == identity['flow_sha256']
surface = pv.read(config['flow']).extract_surface(algorithm='dataset_surface')
surface.points = pose(surface.points)
bounds = tuple(surface.bounds)
corners = pv.Box(bounds=bounds).points
center = np.array(surface.center)
radius = np.linalg.norm(np.ptp(surface.points, axis=0))/2
point_positions = pose(active_x[visible])
ports = {}
for role in ['INLET', 'OUTLET_01', 'OUTLET_02', 'OUTLET_03']:
    path = BASE/'mesh/SV_MESH/mesh-surfaces'/f'{role}.vtp'
    ports[role] = dict(position_mm=pose(np.array(pv.read(path).center)), sha256=digest(path))
q = identity['population_contract']['source']['Q_in_m3_s']
mean_inlet_speed = q / pv.read(BASE/'mesh/SV_MESH/mesh-surfaces/INLET.vtp').area
dt_ms = identity['population_contract']['nominal_dt_s']*1000
font_path = ROOT/'assets/NotoSansCJKsc-Regular.otf'
titlefont = ImageFont.truetype(str(font_path), 30)
smallfont = ImageFont.truetype(str(font_path), 23)

def project(plot, points):
    matrix = plot.camera.GetCompositeProjectionTransformMatrix(plot.renderer.GetTiledAspectRatio(), -1, 1)
    matrix = np.array([[matrix.GetElement(i,j) for j in range(4)] for i in range(4)])
    h = np.column_stack((points, np.ones(len(points)))) @ matrix.T
    ndc = h[:,:2]/h[:,3,None]
    return np.column_stack(((ndc[:,0]+1)*960, (1-ndc[:,1])*540))

def intersects(a, b):
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]

def annotate_ports(plot, pixels, annotations):
    """Choose nearby empty label positions for each new camera angle."""
    image = Image.fromarray(pixels)
    draw = ImageDraw.Draw(image)
    projected_surface = project(plot, surface.points)
    tree = cKDTree(projected_surface)
    occupied = [a['bbox_px'] for a in annotations]
    occupied.append([0,864,384,1080])  # reference-style orientation widget
    labels = {}
    for role, info in ports.items():
        cap = project(plot, np.array([info['position_mm']]))[0]
        text = 'Inlet' if role == 'INLET' else 'Outlet '+str(int(role[-2:]))
        candidates = []
        for distance in [100,140,180,230]:
            for angle in np.arange(0, 360, 30):
                xy = cap + distance*np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))])
                bb = draw.textbbox(tuple(xy), text, font=smallfont, anchor='mm')
                if bb[0] < 24 or bb[2] > 1896 or bb[1] < 203 or bb[3] > 1056:
                    continue
                padded = [bb[0]-7,bb[1]-7,bb[2]+7,bb[3]+7]
                if any(intersects(padded, other) for other in occupied):
                    continue
                in_box = ((projected_surface[:,0] >= padded[0]) & (projected_surface[:,0] <= padded[2])
                          & (projected_surface[:,1] >= padded[1]) & (projected_surface[:,1] <= padded[3]))
                anchor = np.clip(cap, np.array(bb[:2])-6, np.array(bb[2:])+6)
                line_points = anchor+(cap-anchor)*np.linspace(0,.87,24)[:,None]
                proximity = tree.query(line_points)[0]
                penalty = 10*int(in_box.sum()) + .08*int(np.sum(proximity < 4)) + distance*.002
                candidates.append((penalty,tuple(xy),bb,anchor))
        assert candidates, role
        _,xy,bb,anchor = min(candidates, key=lambda entry: entry[0])
        occupied.append([bb[0]-7,bb[1]-7,bb[2]+7,bb[3]+7])
        draw.text(xy,text,font=smallfont,anchor='mm',fill='white')
        draw.line([tuple(anchor),tuple(cap)],fill='#dce2ed',width=2)
        vector = anchor-cap;vector /= np.linalg.norm(vector)
        side = np.array([-vector[1],vector[0]])
        draw.polygon([tuple(cap),tuple(cap+12*vector+4*side),tuple(cap+12*vector-4*side)],fill='#dce2ed')
        labels[role] = dict(text=text,cap_px=cap.tolist(),text_center_px=list(xy),text_bbox_px=list(bb))
    for j,line in enumerate(['Continuous microbubble perfusion',f'Physical time: {TIME:.3f} s',
                            f'dt = {dt_ms:g} ms',f'Inlet mean speed: {mean_inlet_speed:.3f} m/s ({q*60e6:g} mL/min)']):
        draw.text((40,22+40*j),line,font=titlefont,fill='white')
    return image, labels

base_direction = np.array([-.72,.69,.11]);base_direction /= np.linalg.norm(base_direction)
base_azimuth = np.rad2deg(np.arctan2(base_direction[1],base_direction[0]))
base_elevation = np.rad2deg(np.arcsin(base_direction[2]))
view_specs = [(0,base_elevation),(60,25.),(120,40.),(180,base_elevation),(240,25.),(300,40.)]
views = []
for index,(azimuth_offset,elevation) in enumerate(view_specs,1):
    azimuth = base_azimuth+azimuth_offset
    direction = np.array([math.cos(math.radians(elevation))*math.cos(math.radians(azimuth)),
                          math.cos(math.radians(elevation))*math.sin(math.radians(azimuth)),
                          math.sin(math.radians(elevation))])
    plot = pv.Plotter(off_screen=True, window_size=(1920,1080))
    plot.set_background('#000000')
    plot.add_mesh(surface,color='#b2becf',opacity=.16,smooth_shading=True)
    plot.add_mesh(pv.PolyData(point_positions),color='#3399ff',point_size=9.1,
                  render_points_as_spheres=True,lighting=True)
    plot.camera.position = center+direction*radius*4
    plot.camera.focal_point = center
    plot.camera.up = (0,0,1)
    plot.camera.parallel_projection = True
    plot.camera.parallel_scale = radius
    plot.show(auto_close=False)
    controller,orientation = add_spatial_reference_axes(plot,bounds,frame_color='#FFFFFF')
    forward = np.asarray(plot.camera.direction);forward /= np.linalg.norm(forward)
    right = np.cross(forward,[0,0,1]);right /= np.linalg.norm(right)
    up = np.cross(right,forward)
    safe = np.array([[40.,208.],[1880.,1035.]])
    def content_bounds():
        controller.update()
        box = projected_bounds(plot,corners)
        for a in axis_annotations(controller):
            box[0] = np.minimum(box[0],a['bbox_px'][:2]);box[1] = np.maximum(box[1],a['bbox_px'][2:])
        return box
    for _ in range(35):
        box = content_bounds();factor = float(np.max((box[1]-box[0])/(safe[1]-safe[0])))
        plot.camera.parallel_scale *= factor
        box = content_bounds();delta = safe.mean(axis=0)-box.mean(axis=0)
        shift = (-delta[0]*right+delta[1]*up)*(2*plot.camera.parallel_scale/1080)
        plot.camera.position = np.array(plot.camera.position)+shift
        plot.camera.focal_point = np.array(plot.camera.focal_point)+shift
        box = content_bounds()
        if abs(factor-1) < 1e-9:
            break
    assert np.all(box[0] >= safe[0]-1e-4) and np.all(box[1] <= safe[1]+1e-4), box
    plot.render()
    annotations = axis_annotations(controller)
    screenshot, label_positions = annotate_ports(plot,plot.screenshot(return_img=True),annotations)
    name = f'view_{index:02d}_az{azimuth_offset:03d}_el{elevation:04.1f}_t6p0s'
    screenshot.save(OUT/'figures'/f'{name}.png',dpi=(300,300))
    screenshot.convert('RGB').save(OUT/'figures'/f'{name}.pdf',resolution=300)
    gpu = [line for line in plot.ren_win.ReportCapabilities().splitlines()
           if any(word in line for word in ['vendor string','renderer string','version string'])]
    assert any('NVIDIA' in line for line in gpu)
    views.append(dict(name=name,azimuth_offset_deg=azimuth_offset,azimuth_deg=azimuth,
                      elevation_deg=float(elevation),camera_position=list(plot.camera.position),
                      camera_focal_point=list(plot.camera.focal_point),parallel_scale=plot.camera.parallel_scale,
                      content_bounds_px=box.tolist(),vessel_bounds_px=projected_bounds(plot,surface.points).tolist(),
                      port_labels=label_positions,coordinate_annotations=annotations,render_backend=gpu,
                      displayed_particle_ids=shown_ids.tolist(),png_sha256=digest(OUT/'figures'/f'{name}.png')))
    controller.dispose();orientation.dispose();plot.close()
    print('RENDERED',name,len(shown_ids),'particles',flush=True)

contact = Image.new('RGB',(1920,1740),'black')
draw = ImageDraw.Draw(contact)
for i,view in enumerate(views):
    x=(i%2)*960;y=(i//2)*580
    image=Image.open(OUT/'figures'/f"{view['name']}.png").convert('RGB').resize((960,540))
    # Fit overview tiles without cropping; full-resolution originals above.
    contact.paste(image,(x,y))
    draw.text((x+480,y+545),f"View {i+1}: azimuth +{view['azimuth_offset_deg']} deg, elevation {view['elevation_deg']:.1f} deg",
              font=ImageFont.truetype(str(font_path),19),anchor='mt',fill='white')
contact.save(OUT/'CONTACT_SHEET.png',dpi=(300,300))

active_rows = [row for row in rows if row['active_at_6s']]
summary = dict(candidate_id=0,run_name=RUN,physical_time_s=TIME,window_s=[TIME-HALF_WINDOW,TIME+HALF_WINDOW],
               minimum_mean_speed_mm_s=MIN_SPEED_M_S*1000,source_particle_count=len(rows),
               prior_exclusion_count=len(risks|failures),active_at_snapshot=len(active_rows),
               active_prior_exclusion_count=sum(row['prior_solver_failure'] or row['prior_numerical_risk'] for row in active_rows),
               additionally_hidden_slow_count=sum(row['slow_below_2_mm_s'] and not (row['prior_solver_failure'] or row['prior_numerical_risk']) for row in active_rows),
               displayed_count=len(shown_ids),displayed_particle_ids=shown_ids.tolist(),
               partial_window_particle_ids=partial_windows,
               speed_definition='Time average of 3-D speed magnitude over recorded accepted segments intersecting [5.95,6.05] s',
               partial_window_policy='Use available lifetime overlap; never add zero speed after a terminal event or extrapolate',
               maximum_speed_vs_path_length_error_m_s=maximum_speed_check_error,
               maximum_position_vs_held_velocity_error_m=maximum_position_check_error,
               raw_trajectories_unchanged=True,flow_sha256=identity['flow_sha256'],source_sha256=source_hashes,
               style=dict(frame_color='#FFFFFF',background='#000000',bubble_color='#3399ff',vessel_opacity=.16,
                          frame_opacity=.025,grid_opacity=.28,coordinate_units='candidate 0 mm',marker_size_px=9.1),
               renderer_sha256=digest(Path(__file__)),axes_helper_sha256=digest(Path(__file__).parent/'spatial_axes_reference.py'),
               views=views)
(OUT/'data/manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'data/track_hashes.json').write_text(json.dumps(track_hashes,indent=2)+'\n')
cards=''.join(f'''<section><h2>View {i+1}: azimuth +{v['azimuth_offset_deg']}°, elevation {v['elevation_deg']:.1f}°</h2>
<a href="figures/{v['name']}.png"><img src="figures/{v['name']}.png"></a>
<p><a href="figures/{v['name']}.png">PNG</a> · <a href="figures/{v['name']}.pdf">PDF</a></p></section>''' for i,v in enumerate(views))
(OUT/'OPEN_RESULTS.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>Candidate 0 — physical time 6 s</title>
<style>body{{background:#000;color:#dce2ed;font:17px sans-serif;max-width:1600px;margin:30px auto}}img{{width:100%}}a{{color:#3399ff}}section{{margin:30px 0}}</style>
<h1>Candidate 0 — physical time 6.000 s</h1><p>White frame. Same recorded positions and selection in every view. {len(shown_ids)} particles displayed.
Known solver-failure / numerical-risk tracks are hidden, together with particles whose time-averaged speed is below 2 mm/s in the 5.95–6.05 s window.</p>
<p>Display filtering only; raw data are unchanged. <a href="data/display_selection.csv">Selection and speed values</a> · <a href="data/manifest.json">Provenance</a></p>{cards}''')
print(json.dumps({key:summary[key] for key in ['active_at_snapshot','active_prior_exclusion_count','additionally_hidden_slow_count','displayed_count','partial_window_particle_ids']}),flush=True)

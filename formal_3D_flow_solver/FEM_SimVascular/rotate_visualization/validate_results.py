"""Independent checks of the exported FEM visualization and its numerical data."""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_annotations(case, output, name, camera,
                         surface_file='field_diagnostics/data/wall_wss_si.vtp'):
    """Audit exported annotations independently of their layout implementation."""
    plan = json.loads((output/f'{name}_annotations.json').read_text())
    frames = plan['frames']
    assert len(frames) == len(camera) == 432
    assert plan['reference_size_px'] == [1920,1080]
    assert plan['validation']['text_background_rectangles'] == 0
    ids = [r['id'] for r in frames[0]['labels']]
    assert all([r['id'] for r in f['labels']] == ids for f in frames)
    alphas = np.array([[r['opacity'] for r in f['labels']] for f in frames])
    centers = np.array([[r['center_px'] for r in f['labels']] for f in frames])
    assert np.all((alphas >= 0) & (alphas <= 1+1e-12))
    max_alpha_step = float(np.abs(np.roll(alphas,-1,axis=0)-alphas).max())
    max_motion = float(np.linalg.norm(np.roll(centers,-1,axis=0)-centers,axis=2).max())
    assert max_alpha_step < .1 and max_motion < 32.
    for f in frames:
        visible = [r for r in f['labels'] if r['opacity']>1e-8]
        for j,row in enumerate(visible):
            a = row['bbox_px']
            assert 30<a[0]<a[2]<1695 and 105<a[1]<a[3]<980
            for other in visible[j+1:]:
                if row['id'].split(':')[0] == other['id'].split(':')[0]:
                    continue  # Two positions of one crossfading label.
                b = other['bbox_px']
                assert not (a[0]<b[2] and a[2]>b[0] and a[1]<b[3] and a[3]>b[1])
    # Text may have two simultaneous positions during a crossfade. Their total
    # opacity is one, including at the periodic video boundary.
    for text in {r['text'] for r in frames[0]['labels'] if r['kind'] != 'tick'}:
        columns = [j for j,r in enumerate(frames[0]['labels']) if r['text'] == text]
        assert np.allclose(alphas[:,columns].sum(1),1.,atol=1e-12)
    minimum_ticks = {}
    for axis in 'XYZ':
        counts = [sum(r['kind']=='tick' and r['axis']==axis and r['tick_opacity']>.5
                      for r in f['labels']) for f in frames]
        minimum_ticks[axis] = min(counts)
        assert min(counts)>=2,(name,axis,min(counts))
    minimum_clearance, maximum_tip_error = np.inf, 0.
    if any(row['kind']=='port' for row in frames[0]['labels']):
        wall = pv.read(case/surface_file)
        points = np.asarray(wall.points,dtype=float)*1e6
        faces = wall.faces.reshape(-1,4)[:,1:]
        caps = {('Inlet' if role == 'INLET' else 'Outlet '+role[-2:]):
                np.array(pv.read(case/'SV_MESH/mesh-surfaces'/f'{role}.vtp').center)*1e6
                for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']}
        direction = np.array([np.cos(np.deg2rad(6))/np.sqrt(2)]*2+[np.sin(np.deg2rad(6))])
        right = np.cross([0,0,1],direction);right /= np.linalg.norm(right)
        up = np.cross(direction,right)
        for f,cam in zip(frames,camera):
            matrix = np.array(cam['display_matrix'])
            def project(x):
                q = np.asarray(x)@matrix[:3,:3].T+matrix[:3,3]-cam['center_um']
                # Orthographic projection reconstructed from exported camera.
                scale = 1080*.86/(2*cam['parallel_scale_um'])
                return np.column_stack([864+q@right*scale,540-q@up*scale])
            screen = project(points)
            triangles = screen[faces]
            low,high = triangles.min(1),triangles.max(1)
            for row in f['labels']:
                if row['kind'] != 'port' or row['opacity'] < 1e-8:
                    continue
                box = np.array(row['bbox_px'])
                # No triangle bounding box can overlap a port text rectangle;
                # stronger than testing just the cap or the surface vertices.
                overlaps = np.all(high >= box[:2],axis=1)&np.all(low <= box[2:],axis=1)
                assert not overlaps.any(),(f['frame'],row['text'])
                delta = np.maximum(np.maximum(box[:2]-screen,screen-box[2:]),0.)
                clearance = float(np.linalg.norm(delta,axis=1).min())
                assert clearance >= 22.,(f['frame'],row['text'],clearance)
                minimum_clearance = min(minimum_clearance,clearance)
                expected = project([caps[row['text']]])[0]
                tip_error = float(np.linalg.norm(expected-row['arrow_tip_px']))
                assert tip_error < 1e-9
                maximum_tip_error = max(maximum_tip_error,tip_error)
    return dict(frames_checked=432,periodic_seam_checked=True,
                max_opacity_step=max_alpha_step,max_candidate_position_step_px=max_motion,
                no_opaque_text_backgrounds=True,port_text_overlapping_surface_triangles=0,
                visible_text_inside_view=True,unrelated_text_rectangles_overlapping=0,
                min_port_rectangle_to_surface_vertex_distance_px=(minimum_clearance if np.isfinite(minimum_clearance) else None),
                max_arrow_tip_projection_error_px=maximum_tip_error,
                min_visible_tick_marks_per_axis=minimum_ticks,
                all_titles_and_port_names_continuously_present=True,all_pass=True)


def validate(case, output):
    media = json.loads((output/'MEDIA_VALIDATION.json').read_text())
    assert media['all_pass'] and media['z_axis_vertical_all_frames']
    renderer = HERE/'render_visualization.py'
    assert sha(renderer) == media['script_sha256']
    tree = ast.parse(renderer.read_text())
    assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
                   and n.func.attr=='rectangle' for n in ast.walk(tree))
    for relative, digest in media['source_sha256'].items():
        assert sha(case/relative) == digest
    for movie in media['videos']:
        assert movie['decoded_frames'] == movie['distinct_frames'] == 432
        assert movie['size'] == [1920,1080] and movie['fps'] == 24
        assert movie['duration_s'] == 18.
        assert sha(output/movie['file']) == movie['sha256']

    cameras, annotations = {}, {}
    for name in ['full_vessel','branch_vectors']:
        trace = json.loads((output/f'{name}_camera.json').read_text())
        positions = np.array([r['camera_position_um'] for r in trace])
        centers = np.array([r['center_um'] for r in trace])
        axes = np.array([r['axis_unit'] for r in trace])
        up = np.array([r['camera_up_unit'] for r in trace])
        scales = [r['parallel_scale_um'] for r in trace]
        q = positions-centers
        assert len(trace) == 432
        assert np.all(axes == [0,0,1]) and np.all(up == [0,0,1])
        assert np.all(centers == centers[:1]) and np.ptp(scales) == 0
        assert np.ptp(positions[:,2]) < 1e-12
        assert np.ptp(np.linalg.norm(q[:,:2],axis=1)) < 1e-10
        assert np.allclose(np.diff([r['angle_deg'] for r in trace]), 360/432)
        assert max(abs(value) for r in trace for value in r['projected_bounds']) <= .975001
        for direction in q:
            right = np.cross([0.,0.,1.],direction)
            right /= np.linalg.norm(right)
            assert abs(right @ [0.,0.,1.]) < 1e-14
        cameras[name] = dict(frames=432, z_vertical=True, fixed_axis=[0,0,1],
                             constant_height_radius_center_and_zoom=True, geometry_inside_view=True)
        annotations[name] = validate_annotations(case,output,name,trace)

    field = pv.read(case/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    nodes = field.cells.reshape(-1,5)[:,1:]
    vertices = np.asarray(field.points,dtype=float)
    with (output/'data/velocity_glyph_samples.csv').open() as stream:
        samples = list(csv.DictReader(stream))
    velocity_errors, weight_errors = [], []
    for row in samples:
        ids = nodes[int(row['tetra_id'])]
        xyz = vertices[ids]
        point = np.array([float(row[k]) for k in ['x_um','y_um','z_um']])*1e-6
        tail = np.linalg.solve((xyz[1:]-xyz[0]).T,point-xyz[0])
        weights = np.r_[1-tail.sum(),tail]
        assert np.allclose(weights,.25,rtol=0,atol=1e-10)
        vector = weights @ field['Velocity'][ids]
        exported = np.array([float(row[k]) for k in ['vx_m_s','vy_m_s','vz_m_s']])
        assert np.allclose(vector,exported,rtol=1e-10,atol=1e-14)
        assert np.isclose(np.linalg.norm(vector)*1000,float(row['speed_mm_s']),rtol=1e-10)
        velocity_errors.append(float(np.linalg.norm(vector-exported)))
        weight_errors.append(float(np.abs(weights-.25).max()))
    positions = np.array([[float(row[k]) for k in ['x_um','y_um','z_um']] for row in samples])
    distances = np.linalg.norm(positions[:,None,:]-positions[None,:,:],axis=2)
    np.fill_diagonal(distances,np.inf)
    assert distances.min() > 1.25
    result = dict(all_pass=True,case=str(case),output=str(output),camera_checks=cameras,
                  annotation_checks=annotations,
                  vector_count=len(samples),max_velocity_difference_m_s=max(velocity_errors),
                  max_centroid_barycentric_error=max(weight_errors),minimum_vector_spacing_um=float(distances.min()),
                  inputs_and_video_hashes_match=True)
    (output/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case',type=Path,default=HERE/'input_data')
    parser.add_argument('--output',type=Path,default=HERE/'results')
    arguments = parser.parse_args()
    validate(arguments.case.resolve(),arguments.output.resolve())

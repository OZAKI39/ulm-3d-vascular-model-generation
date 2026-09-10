"""Read-only investigation of archived native failures, with new output paths."""
import csv
import json
from pathlib import Path
import re
import time
import numpy as np
from py_scripts.fluid_physics.common import sha256_file, write_json, now
from py_scripts.single_rbc_benchmark.physics import read_off, geometry
from py_scripts.single_rbc_benchmark.analysis import read_csv
from py_scripts.single_rbc_benchmark.quality import mesh_intersections
from .geometry_checks import checked_mesh, membership, independent_intersections, require_same_frame

ROOT = Path(__file__).resolve().parents[2]
LEGACY = ROOT/'data/single_rbc_benchmark/rbc_shear_20260910'
OLD_RUNS = ROOT/'runs/single_rbc_benchmark/rbc_shear_20260910/solver'
TASKS = ('main_mirheo_1','main_mirheo_2','strict_mirheo','diagnostic_mirheo_full')


def read(path):
    return json.loads(Path(path).read_text())


def frame_identity(row):
    return dict(phase=str(row['phase']), global_step=int(row['step']),
                time_star=float(row['time_star']), strain=float(row['strain']))


def audit_task(task, reference, faces, config, output):
    d = OLD_RUNS/task
    spec, execution = read(d/'spec.json'), read(d/'execution.json')
    raw, moments = read_csv(d/'vertices.csv'), read_csv(d/'moments.csv')
    base_geometry = geometry(reference, faces)
    edges = np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]), axis=1), axis=0)
    ref_lengths = np.linalg.norm(reference[edges[:,1]]-reference[edges[:,0]],axis=1)
    keys = list(dict.fromkeys(zip(raw['phase'], raw['step'])))
    first_mismatch = first_intersection = last_normal = normal_before_intersection = None
    first_flags = {}
    metrics, probes = [], []
    raw_count = confirmed_count = uncertain_count = total_count = ambiguous_count = 0
    for frame_index, (phase, step) in enumerate(keys):
        rows = raw[(raw['phase']==phase)&(raw['step']==step)]
        xyz = np.column_stack([rows[k] for k in ('x','y','z')])
        v, g = checked_mesh(rows['vertex'], xyz, faces, config['geometry']['periodic_length'])
        identity = frame_identity(rows[0])
        for k in ('time_star','strain'):
            if len(np.unique(rows[k])) != 1:
                raise ValueError('MIXED_MEMBRANE_FRAME')
        lengths = np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1)
        force = np.linalg.norm(np.column_stack([rows[k] for k in ('fx','fy','fz')]),axis=1)
        legacy_cross = mesh_intersections(v,faces)
        independent = None
        if legacy_cross['nonadjacent_intersections'] or frame_index == len(keys)-1:
            independent = independent_intersections(v,faces,config['repair']['intersection_tolerance'])
            if independent['confirmed_count'] and first_intersection is None:
                normal_before_intersection = last_normal
                first_intersection = dict(**identity, **independent)
        if not legacy_cross['nonadjacent_intersections'] and g['minimum_face_area'] > 1e-10 and g['signed_volume']*base_geometry['signed_volume'] > 0:
            last_normal = dict(**identity, definition='closed ID-consistent nondegenerate mesh, no saved nonadjacent intersections; not a material or continuous-time pass', independent=independent)
        m = moments[(moments['phase']==phase)&(moments['step']==step)]
        if len(m) != 1:
            raise ValueError('MISSING_SAME_TIME_FLUID_MOMENTS')
        q = m[0]
        metric = dict(**identity, area=g['area'], volume=g['volume'], D=g['D'], theta_deg=g['theta_deg'],
                      area_relative_drift=abs(g['area']/base_geometry['area']-1), volume_relative_drift=abs(g['volume']/base_geometry['volume']-1),
                      minimum_face_area=g['minimum_face_area'], minimum_edge=float(lengths.min()),
                      max_reference_edge_ratio=float(np.max(lengths/ref_lengths)),
                      max_wlc_extension_fraction=float(np.max(lengths/ref_lengths)*spec['config']['mirheo_membrane']['x0']),
                      saved_force_max=float(force.max()), saved_force_rms=float(np.sqrt(np.mean(force**2))),
                      saved_force_displacement_scale=float(force.max()*spec['dt']**2/(spec['config']['mirheo_membrane']['mass']*ref_lengths.min())),
                      fluid_N=int(q['N']), inner_N=int(q['inner_N']), temperature=float(q['temperature']), fluid_max_speed=float(q['max_speed']),
                      fluid_speed_times_dt=float(q['max_speed']*spec['dt']), wall_crossings=int(q['wall_crossings']),
                      legacy_intersection_count=legacy_cross['nonadjacent_intersections'],
                      independent_intersection_count=independent['confirmed_count'] if independent else None)
        metrics.append(metric)
        flags = {
            'area_exceeds_original_2pct':metric['area_relative_drift'] > .02,
            'volume_exceeds_original_2pct':metric['volume_relative_drift'] > .02,
            'WLC_above_90pct_max_length':metric['max_wlc_extension_fraction'] > .9,
            'saved_force_displacement_scale_above_0_1':metric['saved_force_displacement_scale'] > .1,
            'fluid_temperature_outside_original_5pct':abs(metric['temperature']-1) > .05,
            'fluid_particle_number_changed':metric['fluid_N'] != int(moments[0]['N'])}
        for name, flagged in flags.items():
            if flagged and name not in first_flags:
                first_flags[name] = metric
        done = int(step)-(spec['prep_steps'] if phase=='shear' else 0)
        p = d/f'fluid_probes_{phase}_{done:08d}.npz'
        if not p.exists():
            raise ValueError('MISSING_SAME_TIME_PROBES')
        with np.load(p) as saved:
            cloud, strain = saved['probes'].copy(), float(saved['strain'])
            require_same_frame(phase,step,str(saved['phase'].item()) if 'phase' in saved else phase,
                               int(saved['global_step']) if 'global_step' in saved else done+(spec['prep_steps'] if phase=='shear' else 0))
        if not np.isclose(strain,identity['strain'],atol=1e-12):
            raise ValueError('PROBE_MEMBRANE_TIME_MISMATCH')
        points = cloud[:,2:5].copy()
        points[:,:2] += config['geometry']['periodic_length']*np.rint((v.mean(axis=0)[:2]-points[:,:2])/config['geometry']['periodic_length'])
        checked = membership(points,cloud[:,0]==1,v,faces,config['repair']['containment_band'])
        n_raw, n_ok, n_band, n_ambiguous = (int(checked[k].sum()) for k in ('raw_mismatch','confirmed_mismatch','uncertain','ambiguous'))
        raw_count += n_raw; confirmed_count += n_ok; uncertain_count += n_band; total_count += len(points); ambiguous_count += n_ambiguous
        row = dict(**identity, probe_file=str(p), tested_point_frames=len(points),
                   raw_winding_mismatches=n_raw, independently_confirmed_mismatches=n_ok,
                   near_surface_uncertain=n_band, ray_disagreement_or_ambiguity=n_ambiguous)
        probes.append(row)
        if n_ok and first_mismatch is None:
            i = int(np.flatnonzero(checked['confirmed_mismatch'])[0])
            tri = int(checked['nearest_triangle'][i])
            first_mismatch = dict(**row, particle_id=int(cloud[i,1]), membership=int(cloud[i,0]),
                                  winding_inside=bool(checked['winding'][i]), ray_inside=bool(checked['winding'][i]),
                                  position=points[i].tolist(), distance_to_membrane=float(checked['distance'][i]),
                                  nearest_triangle=tri, triangle_vertex_ids=faces[tri].tolist(), triangle_coordinates=v[faces[tri]].tolist())
    log_path = d/'shear_00000.log'
    log = log_path.read_text()
    line = next(x for x in log.splitlines() if 'too many triangle collision candidates' in x)
    match = re.search(r'\(coarse\) \((\d+), max (\d+)\)',line)
    last = metrics[-1]
    upper = min(last['global_step']+spec['sample_steps'],spec['prep_steps']+spec['steps'])
    error = dict(log_path=str(log_path), exact_log_line=line, timestamp_as_logged=line.split()[0],
                 time_zone='native local wall-clock; no timezone in the log; archive execution UTC retained separately',
                 coarse_count=int(match[1]), capacity=int(match[2]), exact_native_step=None,
                 completed_steps_exact=None, observed_completed_step_lower_bound=last['global_step'],
                 failing_step_interval=[last['global_step']+1,upper],
                 convention='1-based failing integration attempt; final saved snapshot is only a lower bound',
                 execution_ended_at_utc=execution['ended_at_utc'], exit_code=execution['exit_code'], timeout=execution['timeout'])
    first_kind = 'INSUFFICIENT_TIME_RESOLUTION'
    if first_mismatch and (first_intersection is None or first_mismatch['global_step'] < first_intersection['global_step']):
        first_kind = 'SAMPLED_MEMBERSHIP_DISCREPANCY_PRECEDES_FIRST_CONFIRMED_SELF_INTERSECTION_OR_ERROR'
    elif first_intersection and (first_mismatch is None or first_intersection['global_step'] < first_mismatch['global_step']):
        first_kind = 'CONFIRMED_SAVED_SELF_INTERSECTION_PRECEDES_FIRST_CONFIRMED_PROBE_DISCREPANCY'
    with (output/f'{task}_frame_metrics.csv').open('w') as f:
        w = csv.DictWriter(f,fieldnames=list(metrics[0])); w.writeheader(); w.writerows(metrics)
    return dict(task=task,source_directory=str(d), source_sha256={x:sha256_file(d/x) for x in ('vertices.csv','moments.csv','spec.json','shear_00000.log','execution.json')},
                last_saved_frame=last, planned_strain=spec['steps']*spec['dt']*.02, common_target_reached=False,
                last_geometry_normal_saved_frame=last_normal, first_confirmed_membership_discrepancy=first_mismatch,
                last_normal_saved_frame_before_first_intersection=normal_before_intersection,
                first_confirmed_nonadjacent_intersection=first_intersection, first_screen_events=first_flags,
                sampling=dict(relaxation_time_interval=spec['prep_steps']//10*spec['dt'],shear_time_interval=spec['sample_steps']*spec['dt'],
                              shear_strain_interval=spec['sample_steps']*spec['dt']*.02),
                probe_summary=dict(tested_point_frames=total_count,raw_winding_mismatches=raw_count,
                                   independently_confirmed_mismatches=confirmed_count,near_surface_uncertain=uncertain_count,
                                   ray_ambiguity=ambiguous_count,per_frame=probes,
                                   scope='sampled point-frames, not whole-population crossing counts; no identity comparison across coordinator handoff'),
                error=error, observed_order=first_kind, solver_process_s=execution['elapsed_monotonic_s'],
                unavailable=['membrane one-step velocity and oldPositions','per-step membrane force history','candidate/fine per-triangle distributions and duplicates','whole fluid same-step pre-bounce states','exact failing step'],
                caution='Saved force channel is copied before integration and may include accumulated prior bounce impulse. F*dt²/(m*edge) is a diagnostic scale, not measured displacement. Neither snapshot-average motion nor max_speed*dt is the actual swept distance.')


def initial_force_audit(reference, faces):
    result = {}
    for solver,task in [('HemoCell','candidate1_native_response'),('Mirheo','membrane_probe_candidate')]:
        p = OLD_RUNS/task/'native_response.csv'
        a = read_csv(p)
        a = a[(a['epsilon']==0)&(a['mode']=='shear')]
        forces = np.column_stack([a[k] for k in ('fx','fy','fz')])
        result[solver] = dict(source=str(p),sha256=sha256_file(p),n=len(a),force_rms=float(np.sqrt(np.mean(np.sum(forces**2,axis=1)))),force_max=float(np.linalg.norm(forces,axis=1).max()),net_force=np.sum(forces,axis=0).tolist())
    normals = np.cross(reference[faces[:,1]]-reference[faces[:,0]],reference[faces[:,2]]-reference[faces[:,0]])
    normals /= np.linalg.norm(normals,axis=1)[:,None]
    edge_faces = {}
    for j,f in enumerate(faces):
        for k in range(3):
            edge_faces.setdefault(tuple(sorted((int(f[k]),int(f[(k+1)%3])))),[]).append(j)
    angles = [np.degrees(np.arccos(np.clip(normals[i]@normals[j],-1,1))) for i,j in edge_faces.values()]
    result['reference_unsigned_dihedral_degrees'] = dict(min=float(np.min(angles)),max=float(np.max(angles)),median=float(np.median(angles)),constant_kantor_theta=6.97)
    result['interpretation'] = 'stress_free selects WLC reference lengths/areas; native DihedralKantor uses one constant theta, not the OFF per-edge curvature. Initial total-force audit is independent of fitting the final shear trajectory.'
    return result


def build(config, output):
    start = time.perf_counter()
    reference, faces = read_off(output/'common_reference.off')
    results = []
    for task in TASKS:
        q = audit_task(task,reference,faces,config,output)
        results.append(q)
        write_json(output/f'{task}_timeline.json',q)
        print(json.dumps(dict(task=task,first_probe=q['first_confirmed_membership_discrepancy'],first_intersection=q['first_confirmed_nonadjacent_intersection']),ensure_ascii=False),flush=True)
    timeline = dict(schema_version=2,recorded_at=now(),scope='CPU reanalysis of immutable archived runs; no new solver results',
                    geometry_tolerance=config['repair']['intersection_tolerance'],near_surface_band=config['repair']['containment_band'],
                    tolerance_basis='absolute common length units: intersection 1e-5, membership 1e-4 (twice native bounce position offset 5e-5); selected before inspecting distances',
                    vertex_check='CSV vertex field is a sorted row index. Archived worker sorts native IDs and checks their unique contiguity when reading same-frame saved forces. Initial raw native IDs and crash-state IDs were not persisted. CSV connectivity, closed orientation and integer-period unwrap pass; no rotation fitting.',
                    runs=results,initial_force_audit=initial_force_audit(reference,faces),cpu_analysis_s=time.perf_counter()-start)
    write_json(output/'failure_timeline.json',timeline)
    return timeline

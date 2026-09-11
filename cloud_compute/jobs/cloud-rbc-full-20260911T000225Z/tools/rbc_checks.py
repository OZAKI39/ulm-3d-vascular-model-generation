"""CPU-only saved-frame screening shared by the live guard and final report."""
from common import *
import csv, math, re
import numpy as np
from cloud_geometry import read_off, geometry
from rbc_geometry import checked_mesh, independent_intersections, membership
from rbc_h5 import read_xmf


def frame_check(ids, positions, velocities, faces, reference, spec):
    c = spec['config']; repair = c['repair']
    v, g = checked_mesh(ids, positions, faces, c['geometry']['periodic_length'])
    ref = geometry(reference, faces)
    edges = np.unique(np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0)
    ratios = np.linalg.norm(v[edges[:, 1]] - v[edges[:, 0]], axis=1) / np.linalg.norm(reference[edges[:, 1]] - reference[edges[:, 0]], axis=1)
    extension = float(ratios.max() * c['mirheo_membrane']['x0'])
    intersections = independent_intersections(v, faces, repair['intersection_tolerance'])
    hard = []
    if not np.isfinite(velocities).all(): hard.append('NONFINITE_MEMBRANE_VELOCITY')
    if g['minimum_face_area'] <= 1e-12: hard.append('DEGENERATE_TRIANGLE')
    if extension >= 1: hard.append('WLC_EXTENSION_OUTSIDE_DOMAIN')
    if intersections['confirmed_count']: hard.append('CONFIRMED_SELF_INTERSECTION')
    if v[:, 2].min() < -1e-5 or v[:, 2].max() > c['geometry']['gap'] + 1e-5: hard.append('MEMBRANE_WALL_CROSSING')
    area = abs(g['area']/ref['area'] - 1); volume = abs(g['volume']/ref['volume'] - 1)
    return dict(geometry=g, area_relative_drift=area, volume_relative_drift=volume,
                max_WLC_extension=extension, max_membrane_speed=float(np.linalg.norm(velocities, axis=1).max()) if np.isfinite(velocities).all() else None,
                intersections=intersections, hard_failures=hard,
                soft_failures=[name for name, val in [('area_relative_drift', area), ('volume_relative_drift', volume)] if val > c['screen'][name]],
                vertices=v.tolist(), vertex_ids=np.sort(ids).tolist())


def scan_frames(root, spec, cache=None, final=False):
    """A native indexed frame is a progress lower bound, never a run return."""
    root = Path(root); cache = {} if cache is None else cache
    ref, faces = read_off(spec['mesh']); pad = spec['config']['geometry']['wall_layer']
    for phase, interval in [('relaxation', spec['prep_sample_steps']), ('shear', spec['sample_steps'])]:
        paths=sorted((root/f'native_{phase}').glob('rbc*.xmf'))
        returned=(root/f'phase_{phase}.json').exists() and read(root/f'phase_{phase}.json').get('native_returned')
        if not final and not returned: paths=paths[:-1]  # Next native frame proves the previous writer has advanced/closed.
        for path in paths:
            key = path.relative_to(root).as_posix()
            if key in cache: continue
            step = int(re.search(r'(\d+)\.xmf$', path.name)[1]) * interval
            try:
                raw = read_xmf(path)
            except Exception as ex:
                # A file may still be in flight while the solver is alive.
                if final: cache[key] = dict(phase=phase, phase_step=step, unreadable=str(ex), hard_failures=[])
                continue
            try:
                check = frame_check(raw['id'], raw['position']-[0,0,pad], raw['velocity'], faces, ref, spec)
            except Exception as ex:
                check = dict(hard_failures=['INVALID_GEOMETRY:'+str(ex)])
            check.update(phase=phase, phase_step=step, step=step+(spec['prep_steps'] if phase=='shear' else 0),
                         time_star=step*spec['dt'], strain=step*spec['dt']*spec['config']['protocol']['shear_rate'] if phase=='shear' else 0,
                         source=key, observation='native beforeForces; confirmed lower bound')
            cache[key] = check
        for path in sorted(root.glob(f'membrane_state_{phase}_*.npz')):
            key = path.name
            if key in cache: continue
            try:
                with np.load(path, allow_pickle=False) as raw:
                    step = int(path.stem.rsplit('_',1)[1])
                    check = frame_check(raw['ids'],raw['positions'],raw['velocities'],faces,ref,spec)
                check.update(phase=phase,phase_step=step,step=step+(spec['prep_steps'] if phase=='shear' else 0),time_star=step*spec['dt'],
                             strain=step*spec['dt']*spec['config']['protocol']['shear_rate'] if phase=='shear' else 0,
                             source=key,observation='after successful native return')
                cache[key] = check
            except Exception as ex:
                if final: cache[key] = dict(phase=phase,phase_step=step,unreadable=str(ex),hard_failures=[])
    return cache


def preparation(frames, reference, spec):
    """Same nonoverlapping windows/residual/axial-angle definitions as preparation_quality.py."""
    criteria = spec['preparation_criteria']; width = criteria['window_time_star']
    samples = sorted([x for x in frames if x.get('phase')=='relaxation' and 'vertices' in x], key=lambda x:x['phase_step'])
    ref = reference-reference.mean(0); scale = spec['reference_radius']; windows=[]; transitions=[]
    for n in range(int(spec['prep_steps']*spec['dt']/width+1e-9)):
        lo=n*width;hi=lo+width; selected=[f for f in samples if lo-1e-9 <= f['time_star'] < hi-1e-9]
        if not selected: continue
        coords=[np.array(f['vertices']) for f in selected]; angles=[f['geometry']['theta_deg'] for f in selected]
        mean=np.mean([v-v.mean(0) for v in coords],axis=0)
        angle=None if any(a is None for a in angles) else float(np.degrees(.5*np.angle(np.mean(np.exp(2j*np.radians(angles))))))
        windows.append(dict(begin=lo,end=hi,samples=len(selected),coordinates=mean,angle_deg=angle,
                            reference_shape_rms_over_a=float(np.sqrt(np.mean(np.sum((mean-ref)**2,axis=1)))/scale)))
    for before, after in zip(windows, windows[1:]):
        residual=float(np.sqrt(np.mean(np.sum((after['coordinates']-before['coordinates'])**2,axis=1)))/scale/width)
        angle=None if before['angle_deg'] is None or after['angle_deg'] is None else abs((after['angle_deg']-before['angle_deg']+90)%180-90)/width
        transitions.append(dict(end=after['end'],residual=residual,angle_change_per_time=angle,
                                passed=residual<=criteria['residual_shape_rms_over_a_per_time'] and angle is not None and angle<=criteria['orientation_change_degrees_per_time']))
    required=criteria['minimum_consecutive_windows']; last=transitions[-required:]
    enough=bool(samples) and samples[-1]['phase_step']==spec['prep_steps'] and len(windows)==int(spec['prep_steps']*spec['dt']/width+1e-9)
    residual_pass=enough and len(last)==required and all(t['passed'] for t in last)
    geometry_pass=bool(samples) and all(not f['hard_failures'] and not f['soft_failures'] for f in samples)
    for w in windows: w.pop('coordinates')
    return dict(passed=residual_pass and geometry_pass,status='PASS' if residual_pass and geometry_pass else 'PREPARATION_QUALITY_FAILED_OR_INCOMPLETE',
                criteria=criteria,complete_window_evidence=enough,residual_windows_pass=residual_pass,geometry_pass=geometry_pass,windows=windows,transitions=transitions,
                equilibrium_certified=False,reference_shape_comparison='diagnostic only; no new HemoCell matching test',decision='shear only after all preparation gates pass; no extension past t*=30')


def probe_checks(root, spec, frames):
    root=Path(root);_,faces=read_off(spec['mesh']);L=spec['config']['geometry']['periodic_length'];pad=spec['config']['geometry']['wall_layer'];rows=[]
    lookup={(f['phase'],f['phase_step']):f for f in frames if 'vertices' in f}
    inputs=[]
    for p in sorted(root.glob('fluid_probes_*.npz')):
        phase=p.stem.split('_')[2];step=int(p.stem.rsplit('_',1)[1])
        with np.load(p,allow_pickle=False) as raw: inputs.append((phase,step,raw['probes'].copy(),p.name,'96 evenly spaced sorted IDs per side, chosen at returned endpoint'))
    for phase in ['relaxation','shear']:
        paths=[root/f'native_{phase}'/f'{name}00000.xmf' for name in ['outer','inner']]
        if not all(p.exists() for p in paths): continue
        try:
            probes=[]
            for member,p in enumerate(paths):
                r=read_xmf(p);order=np.argsort(r['id']);take=order[np.arange(0,len(order),max(1,len(order)//96))[:96]]
                if not np.isfinite(r['position']).all() or not np.isfinite(r['velocity']).all():raise ValueError('NONFINITE_INITIAL_FLUID')
                probes.extend([[member,int(r['id'][i]),*(r['position'][i]-[0,0,pad])] for i in take])
            inputs.append((phase,0,np.array(probes),','.join(p.name for p in paths),'96 evenly spaced sorted IDs per side at native initial frame; no tracked path between frames'))
        except Exception as ex: rows.append(dict(phase=phase,phase_step=0,status='UNREADABLE',error=str(ex)))
    for phase,step,probes,source,rule in inputs:
        frame=lookup.get((phase,step))
        if frame is None: rows.append(dict(phase=phase,phase_step=step,status='MEMBRANE_FRAME_MISSING'));continue
        vertices=np.array(frame['vertices']);points=probes[:,2:5].copy();points[:,:2]+=L*np.rint((vertices.mean(0)[:2]-points[:,:2])/L)
        try:
            m=membership(points,probes[:,0]==1,vertices,faces,spec['config']['repair']['containment_band'])
            rows.append(dict(phase=phase,phase_step=step,status='CHECKED',tested=len(points),selection=rule,source=source,
                             raw_mismatches=int(m['raw_mismatch'].sum()),confirmed_mismatches=int(m['confirmed_mismatch'].sum()),
                             near_surface_uncertain=int(m['uncertain'].sum()),second_method_ambiguous=int(m['ambiguous'].sum()),
                             confirmed_mismatch_ids=probes[m['confirmed_mismatch'],1].astype(int).tolist()))
        except Exception as ex: rows.append(dict(phase=phase,phase_step=step,status='UNREADABLE',error=str(ex)))
    return dict(strict_impermeability='NOT_VERIFIED',scope='At most 192 probes per saved initial/returned phase endpoint; no claim about unsaved times or all fluid particles.',
                checks=rows,tested_point_frames=sum(r.get('tested',0) for r in rows),confirmed_mismatches=sum(r.get('confirmed_mismatches',0) for r in rows),
                near_surface_uncertain=sum(r.get('near_surface_uncertain',0) for r in rows),classifier_corrections_within_phase=0,
                handoff_classifier='native initial classification after sorted-old-ID remapping; disclosed separately')


def outcome(spec, completion, exit_code, frames):
    reached=bool(completion and completion.get('completed') and completion.get('shear_steps')==spec['steps'] and
                 completion.get('prep_steps')==spec['prep_steps'] and exit_code==0 and abs(completion.get('strain_end',-1)-4)<1e-5 and
                 any(f.get('phase')=='shear' and f.get('phase_step')==spec['steps'] and 'vertices' in f for f in frames))
    progress={phase:max([f.get('phase_step',0) for f in frames if f.get('phase')==phase],default=0) for phase in ['relaxation','shear']}
    return dict(CLOUD_RBC_RUN_COMPLETE='PASS' if reached else 'NOT_COMPLETE',actual_successful_returned_prep_steps=completion.get('prep_steps',0) if completion else None,
                actual_successful_returned_shear_steps=completion.get('shear_steps',0) if completion else None,
                last_saved_progress_lower_bound=progress,strain_lower_bound=progress['shear']*spec['dt']*spec['config']['protocol']['shear_rate'],
                PHYSICAL_MODEL_VALIDATION='NOT_MATCHED; material calibration, space/time convergence and HemoCell comparison NOT_TESTED this round',qualified_speedup=None)

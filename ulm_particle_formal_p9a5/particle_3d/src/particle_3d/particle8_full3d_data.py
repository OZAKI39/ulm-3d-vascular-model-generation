"""Read-only input discovery and camera-independent trajectory extraction for P8 3D."""
from pathlib import Path
import subprocess
import math
import numpy as np
from .particle8_replay import REPO,SCHEMA,read,write,digest,canonical_hash,snapshot,interpolate,validate_scene
from .audit import read_frozen

OUTPUT=REPO/'particle_3d/outputs/particle8_full3d'
UPSTREAM=REPO/'particle_3d/reports/particle8'
CASES={
 'A_real_single_mb_orbit':dict(scene='real_single_mb_entry',camera='orbit',classification='REAL FROZEN GEOMETRY / SINGLE MB / SAVED P6.5 VIA P7 / NO RBC',title='Single MB entry in the full Frozen vessel',tail_window_s=.00035),
 'A_real_single_mb_fixed':dict(scene='real_single_mb_entry',camera='fixed',classification='REAL FROZEN GEOMETRY / SINGLE MB / SAVED P6.5 VIA P7 / NO RBC',title='Single MB entry | fixed-camera reference',tail_window_s=.00035),
 'B_real_mixed_smoke_orbit':dict(scene='real_mixed_inlet_smoke',camera='orbit',classification='REAL FROZEN GEOMETRY / SMOKE',title='Admitted particles in the full Frozen vessel',tail_window_s=.00035),
 'C_synthetic_lifecycle_orbit':dict(scene='synthetic_open_section',camera='orbit',classification='SYNTHETIC CONTROL / NOT REAL FROZEN LUMEN',title='3D lifecycle and outlet deletion',tail_window_s=.00025),
 'D_restart_continuity_orbit':dict(scene='synthetic_open_section',camera='orbit',classification='SYNTHETIC CONTROL / RESTART AUDIT',title='Checkpoint continuity in a rotating 3D view',tail_window_s=.00025),
}
LABELS=[
 'DISPLAY-ONLY INTERPOLATION between saved states; no new solver states.',
 'Camera motion is display-only; physical coordinates and event times are unchanged.',
 'Real RBC passage NOT established | H_D = 0.45 is feed volume fraction; tube Hct is a separate diagnostic.',
 'Production timestep / neighbor settings NOT FROZEN | sphere-normal near-field only; non-spherical lubrication NOT FROZEN.',
 'No Particle-7.5 | No CFD | No new RBC deformation physics | No full suspension / PK.',
]


def safe_output(path):
    out=Path(path).resolve()
    old=REPO/'particle_3d/reports'
    if out==old or old in out.parents or out in old.parents:
        raise ValueError('New outputs must not overlap existing scientific reports')
    return out


def old_files():
    """Lock every tracked pre-existing particle file, not just selected screenshots."""
    paths=subprocess.check_output(['git','ls-files','particle_3d'],cwd=REPO,text=True).splitlines()
    return [p for p in paths if '/outputs/particle8_full3d/' not in p and 'full3d' not in p.lower()]


def verify_preservation(root):
    lock=read(Path(root)/'data/upstream_lock.json')
    for path,h in lock['sha256'].items():
        if digest(REPO/path)!=h:raise ValueError('Upstream file changed: '+path)
    return len(lock['sha256'])


def discover():
    scenes={}
    for p in sorted((REPO/'particle_3d/reports').glob('particle8/**/**scene*.json')):
        d=read(p)
        if d.get('schema')==SCHEMA:
            validate_scene(d)
            if d['name'] in scenes and scenes[d['name']]['sha256']!=digest(p):raise ValueError('Conflicting upstream scenes')
            scenes[d['name']]=dict(path=str(p.relative_to(REPO)),sha256=digest(p),source_classification=d['source_classification'],
                records=len(d['records']),time_range_s=d['time_range_s'])
    for name in set(c['scene'] for c in CASES.values())|{'synthetic_restarted'}:
        if name not in scenes:raise ValueError('Required replay scene not discovered: '+name)
    return scenes


def prepare(output=OUTPUT):
    root=safe_output(output);(root/'data').mkdir(parents=True,exist_ok=True)
    lock=root/'data/upstream_lock.json'
    if lock.exists():verify_preservation(root)
    else:write(lock,dict(baseline_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
                        sha256={p:digest(REPO/p) for p in old_files()}))
    scenes=discover();write(root/'data/discovered_inputs.json',dict(scenes=scenes,
        scripts=sorted(str(p.relative_to(REPO)) for p in (REPO/'particle_3d/scripts').glob('*particle8*.py') if 'full3d' not in p.name),
        upstream_output=str(UPSTREAM.relative_to(REPO)),new_output=str(root)))
    provenance,mesh,flow,boundaries=read_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    arrays={};summary={}
    for role,surf in boundaries.items():
        faces=np.asarray(surf.faces).reshape(-1,4)
        assert np.all(faces[:,0]==3)
        arrays[role+'_points_m']=np.asarray(surf.points,dtype=float)
        arrays[role+'_faces']=faces[:,1:].astype(np.int64)
        summary[role]=dict(points=surf.n_points,triangles=surf.n_cells,
            bounds_m=list(map(float,surf.bounds)),source=provenance['boundary_manifest']['boundaries'][role])
    np.savez_compressed(root/'data/full_frozen_geometry.npz',**arrays)
    write(root/'data/full_geometry_manifest.json',dict(source='OFFICIAL_FROZEN_BOUNDARY_MANIFEST',
        surfaces=summary,full_wall_triangle_count=boundaries['WALL'].n_cells,decimation=False,physical_coordinate_changes=False,
        source_mesh_sha256=provenance['mesh_sha256'],source_flow_sha256=provenance['flow_sha256'],units='m'))
    write(root/'data/display_contract.json',dict(cases=CASES,labels=LABELS,full3d=True,projection='ORTHOGRAPHIC',
        camera=dict(start_azimuth_deg=35,end_azimuth_deg=115,elevation_deg=28,real_inset_azimuth_offset_deg=90,
            role='DISPLAY_ONLY_CAMERA_ORBIT_NO_PHYSICAL_COORDINATE_ROTATION'),
        synthetic_display=dict(z_exaggeration=20000,physical_length_m=1e-8,physical_width_m=200e-6,
            main_view='CENTER_GLYPHS_NOT_FINITE_SHAPE_GEOMETRY',secondary_view='C_UNSTRETCHED_ORIGINAL_FINITE_SHAPES; D_RESTARTED_CENTER_GLYPHS_SAME_STRETCH'),
        real_display=dict(unit_conversion='(position_m - fixed_origin_m) * 1e6',axis_stretching=False,
            main_view='ALL_OFFICIAL_WALL_FACES_AND_ALL_FOUR_OPEN_BOUNDARIES',secondary_view='LOCAL_ZOOM_SAME_PHYSICAL_SHAPE'),
        tail_mode_default='recent',tail_geometry='SAVED_POLYLINE_CLIPPED_AT_RECENT_TIME_AND_CURRENT_TIME_WITH_LINEAR_DISPLAY_ENDPOINTS',
        tail_visibility='ACTIVE_PARTICLES_ONLY_REMOVE_WITH_PARTICLE_AT_DELETE',pending_display='COUNTS_ONLY_NO_SPATIAL_POSITIONS'))
    return scenes


def scene_from_name(root,name):
    spec=read(Path(root)/'data/discovered_inputs.json')['scenes'][name]
    path=REPO/spec['path']
    if digest(path)!=spec['sha256']:raise ValueError('Scene source changed')
    return read(path)


def trajectory_tail(record,t,window_s=None):
    """Only original polyline knots plus clipped endpoints; no spline or future tail."""
    if record['admit_time_s'] is None or t<record['admit_time_s'] or (record['delete_time_s'] is not None and t>=record['delete_time_s']):return []
    if window_s is not None and (not math.isfinite(window_s) or window_s<=0):raise ValueError('Positive finite tail window required')
    start=max(record['admit_time_s'],t-window_s) if window_s is not None else record['admit_time_s']
    times=[start]+[s['time_s'] for s in record['trajectory'] if start<s['time_s']<t]
    if t>start:times.append(t)
    return [dict(time_s=float(u),**interpolate(record,float(u))) for u in times]


def schedule(case,scene):
    if case.startswith('A_'):return [(float(t),'ENTRY REPLAY') for t in np.linspace(*scene['time_range_s'],180)]
    if case.startswith('B_'):
        end=scene['time_range_s'][1]
        return [(float(t),'ADMISSION ONLY / STATIC OBSTACLE') for t in np.r_[np.linspace(0,end,165),np.repeat(end,15)]]
    if case.startswith('C_'):
        windows=[(0,.002,'STARTUP'),(.117,.1184,'CUT TO FIRST MB WINDOW'),(.2347,.2361,'CUT TO SECOND MB WINDOW')]
        return [(float(t),label) for lo,hi,label in windows for t in np.linspace(lo,hi,60)]+[(.25,'CUT TO FINAL ACCOUNTING')]*15
    return [(float(t),'CHECKPOINT AT 0.125 s') for t in np.r_[np.linspace(.124,.125,70,endpoint=False),np.repeat(.125,10),np.linspace(.125,.126,70)]]+[(.25,'CUT TO FINAL ACCOUNTING')]*10


def camera_angles(case,frame,total):
    phi=35. if CASES[case]['camera']=='fixed' else 35.+80.*frame/max(1,total-1)
    return dict(azimuth_deg=phi,elevation_deg=28.,role='DISPLAY_ONLY')


def replay_frame(scene,t,tail_window_s,tail_mode='recent'):
    if tail_mode not in ['recent','full']:raise ValueError('Unknown tail mode')
    raw=snapshot(scene,t);records={r['particle_id']:r for r in scene['records']}
    particles=[]
    for p in raw['active']:
        r=records[p['particle_id']]
        particles.append(dict(**p,geometry_sha256=canonical_hash(r['geometry']),admitted_shape_sha256=canonical_hash(r['admitted_shape']),
            tail=trajectory_tail(r,t,None if tail_mode=='full' else tail_window_s)))
    # Pending identity remains in the upstream scene; never allocate a render actor.
    return dict(time_s=t,source_classification=scene['source_classification'],counts=raw['counts'],
        restart_segment=raw['restart_segment'],particles=particles,
        pending_ids=[p['particle_id'] for p in raw['pending']],pending_drawn_ids=[],
        tail_mode=tail_mode,tail_window_s=tail_window_s if tail_mode=='recent' else None)

"""Shared saved-trajectory scene for P8.1 figures, animations and audits."""
from pathlib import Path
from collections import Counter
import csv,sys,platform
from importlib.metadata import version
import numpy as np
from .particle8_replay import read,digest,canonical_hash,slerp
from .particle81_simulation import OUTPUT,STATES,dump,lock_upstream

TITLE='Frozen flow full-vessel ULM-style microbubble trajectory simulation'
INDEPENDENT='ULM replay is assembled from independently integrated full-vessel microbubble trajectories in a frozen flow field.'
FOOTERS=[
 'Complete MB trajectory visualization layer | Independent superposition | No RBC-coupled physiological suspension',
 'Real RBC passage unresolved | No full suspension / PK | No new CFD physics / P7.5 / RBC deformation',
 'H_D = 0.45 feed fraction; tube Hct separate | C_MB = 8.5e12 m^-3 nominal | isotropic orientation V0 assumption',
 'P6.5 sphere-normal regularization only | production dt / neighbors / nonspherical lubrication NOT FROZEN',
]
OUTLETS=['OUTLET_01','OUTLET_02','OUTLET_03']


def assemble(root=OUTPUT):
    root=Path(root);births=read(root/'data/birth_ledger.json');entries=[];events=[];checks=[]
    fonts=[Path('/usr/share/fonts/truetype/dejavu')/name for name in ['DejaVuSans.ttf','DejaVuSans-Bold.ttf']]
    dump(root/'data/software_environment.json',dict(python=sys.version,platform=platform.platform(),
        packages={name:version(name) for name in ['numpy','scipy','pyvista','vtk','matplotlib','Pillow','imageio-ffmpeg']},
        font_sha256={str(p):digest(p) for p in fonts},
        pixel_determinism_scope='REPEATED_RENDER_IN_THIS_SAVED_SOFTWARE_ENVIRONMENT; NOT_A_CROSS_GPU_GUARANTEE'))
    for birth in births['events']:
        pid=birth['particle_id'];meta_path=root/'trajectories'/f'mb_{pid:06d}.json';meta=read(meta_path)
        path=root/meta['samples_path'];data=np.load(path)['samples']
        if digest(path)!=meta['samples_sha256']:raise ValueError('Saved trajectory hash mismatch')
        if canonical_hash(birth)!=meta['birth_metadata_sha256']:raise ValueError('Birth metadata mismatch')
        dt=np.diff(data[:,0]) if len(data) else np.array([])
        check=dict(particle_id=pid,finite=bool(np.isfinite(data).all()),strictly_monotone=bool(np.all(dt>0)),
            sample_count=len(data),completed=meta['completed'],end_reason=meta['end_reason'],
            original_q_preserved=bool(not len(data) or np.array_equal(data[0,10:14],birth['q'])),
            displacement_equation_max_error_m=float(np.max(np.abs(np.diff(data[:,1:4],axis=0)-dt[:,None]*data[1:,4:7]))) if len(data)>1 else 0.,
            minimum_g_nf_m=float(data[:,15].min()) if len(data) else None)
        check['pass_checks']=check['finite'] and check['strictly_monotone'] and check['original_q_preserved']
        if meta['completed']:
            check['pass_checks'] &= len(data)>=2 and meta['exit_outlet'] in OUTLETS and meta['birth_time_s']<=meta['exit_time_s']
        checks.append(check)
        entry={key:meta.get(key) for key in ['particle_id','birth_time_s','admit_time_s','radius_m','diameter_um','initial_q',
            'inlet_face','inlet_barycentric','inlet_position_m','completed','exit_outlet','exit_time_s','residence_time_s',
            'path_length_m','end_reason','failure_detail','sample_count','last_elapsed_time_s','last_physical_time_s',
            'minimum_original_wall_gap_m','maximum_refinement_depth','rejected_trials','nearfield_states','integration_config']}
        entry.update(samples_path=meta['samples_path'],samples_sha256=meta['samples_sha256'],
            metadata_path=str(meta_path.relative_to(root)),metadata_sha256=digest(meta_path))
        entries.append(entry)
        events.append(dict(particle_id=pid,time_s=birth['birth_time_s'],kind='SCHEDULED_BIRTH',order=0))
        if len(data):
            events.append(dict(particle_id=pid,time_s=birth['birth_time_s'],kind='ADMITTED_ACTIVE',order=1))
            events.append(dict(particle_id=pid,time_s=meta['last_physical_time_s'],kind='OUTLET_EXIT' if meta['completed'] else 'EXPLICIT_SAFETY_STOP',
                reason=meta['end_reason'],outlet=meta['exit_outlet'],order=2))
            if meta['completed']:events.append(dict(particle_id=pid,time_s=meta['exit_time_s'],kind='DELETED_AFTER_EXIT',order=3))
        else:events.append(dict(particle_id=pid,time_s=birth['birth_time_s'],kind='UNRESOLVED_ADMISSION_TERMINAL_RECORD',reason=meta['end_reason'],order=2))
    events.sort(key=lambda e:(e['time_s'],e['particle_id'],e['order']))
    completed=[e for e in entries if e['completed']];outlet_counts={o:sum(e['exit_outlet']==o for e in completed) for o in OUTLETS}
    reps=[]
    for outlet in OUTLETS:
        cohort=[e for e in completed if e['exit_outlet']==outlet]
        if cohort:
            cohort.sort(key=lambda e:(e['residence_time_s'],e['particle_id']));reps.append(cohort[len(cohort)//2]['particle_id'])
    for e in sorted(completed,key=lambda e:e['path_length_m'],reverse=True):
        if len(reps)>=3:break
        if e['particle_id'] not in reps:reps.append(e['particle_id'])
    cat=dict(schema='PARTICLE81_SAVED_FULL_TRAJECTORY_SCENE_V1',title=TITLE,classification='REAL_FROZEN_VESSEL_AND_FLOW_INDEPENDENT_MB',
        independent_superposition=INDEPENDENT,units=dict(time='s',position='m',velocity='m/s',omega='1/s',
            radius='m',diameter_um='um',gap='m',path_length='m',quaternion='dimensionless'),original_sample_columns=read(root/entries[0]['metadata_path'])['columns'],
        state_codes=STATES,entries=entries,representative_ids=reps,selection_rule='Median residence in each observed outlet; fill with longest distinct completed paths',
        acquisition_birth_window_s=births['acquisition_birth_window_s'],replay_end_time_s=max(e.get('last_physical_time_s') or e['birth_time_s'] for e in entries),
        scheduled=len(entries),admitted=sum(e['sample_count']>0 for e in entries),completed=len(completed),
        end_reasons=dict(Counter(e['end_reason'] for e in entries)),outlet_counts=outlet_counts,
        completed_path_length_range_m=[min(e['path_length_m'] for e in completed),max(e['path_length_m'] for e in completed)] if completed else [],
        physical_samples=sum(e['sample_count'] for e in entries),scientific_labels=FOOTERS,
        source_hashes=dict(birth_ledger=digest(root/'data/birth_ledger.json'),frozen_geometry=digest(root/'data/full_frozen_geometry.npz'),
            upstream_lock=digest(root/'data/upstream_lock.json')),
        display_only=dict(position_interpolation='LINEAR_BETWEEN_ACCEPTED_SUBSTEPS',quaternion_interpolation='P8_SLERP',
            extrapolation=False,physical_axis_stretching=False,renderer_no_rng=True,paths='SAVED_SOLVER_POLYLINES_NO_SPLINE',
            history='COMPLETED_TRACKS_ONLY_IN_ACCUMULATION; ACTIVE_PARTIAL_TAILS_IN_REPLAY',
            active_counts='WITHIN_SAVED_COMPUTATION_INTERVALS_ONLY',stopped_particle_fate='UNKNOWN_AFTER_COMPUTATIONAL_CENSORING'))
    dump(root/'data/particle8_1_trajectory_catalog.json',cat)
    dump(root/'data/particle8_1_events.json',dict(events=events,clock='ABSOLUTE_ACQUISITION_SECONDS',tie_order='BIRTH_ADMIT_EXIT_DELETE'))
    dump(root/'data/trajectory_audit.json',dict(all_pass=all(c['pass_checks'] for c in checks),trajectories=checks,
        minimum_completed_required=300,preferred_completed=1000,minimum_scale_reached=len(completed)>=300,
        preferred_scale_reached=len(completed)>=1000,all_three_outlets_observed=all(outlet_counts.values()),
        upstream_preserved_files=lock_upstream()))
    with (root/'data/trajectory_catalog.csv').open('w',newline='') as f:
        keys=['particle_id','birth_time_s','diameter_um','completed','exit_outlet','exit_time_s','residence_time_s','path_length_m','sample_count','end_reason']
        writer=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');writer.writeheader();writer.writerows({k:e[k] for k in keys} for e in entries)
    return cat


class Scene:
    def __init__(self,root=OUTPUT,verify=True):
        self.root=Path(root);self.catalog=read(self.root/'data/particle8_1_trajectory_catalog.json')
        self.sha256=digest(self.root/'data/particle8_1_trajectory_catalog.json')
        if verify:
            for key,name in [('birth_ledger','birth_ledger.json'),('frozen_geometry','full_frozen_geometry.npz'),('upstream_lock','upstream_lock.json')]:
                if digest(self.root/'data'/name)!=self.catalog['source_hashes'][key]:raise ValueError('Scene source changed: '+name)
        self.entries={e['particle_id']:e for e in self.catalog['entries']};self.arrays={}
        for pid,e in self.entries.items():
            path=self.root/e['samples_path']
            if verify and digest(path)!=e['samples_sha256']:raise ValueError('Trajectory source changed')
            if verify and digest(self.root/e['metadata_path'])!=e['metadata_sha256']:raise ValueError('Trajectory metadata changed')
            data=np.load(path)['samples'];data.flags.writeable=False;self.arrays[pid]=data
        self.completed=[pid for pid,e in self.entries.items() if e['completed']]

    def position(self,pid,elapsed):
        a=self.arrays[pid]
        if len(a)<1 or elapsed<a[0,0] or elapsed>a[-1,0]:raise ValueError('No replay extrapolation')
        k=int(np.searchsorted(a[:,0],elapsed,side='right'))-1
        if k>=len(a)-1 or a[k,0]==elapsed:return a[k,1:4].copy(),a[k,10:14].copy()
        weight=(elapsed-a[k,0])/(a[k+1,0]-a[k,0])
        return (1-weight)*a[k,1:4]+weight*a[k+1,1:4],np.array(slerp(a[k,10:14],a[k+1,10:14],weight))

    def snapshot(self,t,only_ids=None):
        active=[];finished=[];counts=Counter(scheduled=0,admitted=0,active=0,completed=0,unresolved=0,stopped=0)
        for pid,e in self.entries.items():
            if only_ids is not None and pid not in only_ids:continue
            if t<e['birth_time_s']:continue
            counts['scheduled']+=1
            if e['sample_count']==0:counts['unresolved']+=1;continue
            counts['admitted']+=1
            if t>=e['last_physical_time_s']:
                if e['completed']:counts['completed']+=1;finished.append(pid)
                else:counts['stopped']+=1
            else:
                tau=float(np.clip(t-e['birth_time_s'],0,self.arrays[pid][-1,0]));position,q=self.position(pid,tau)
                active.append(dict(particle_id=pid,elapsed_time_s=tau,position_m=position.tolist(),q=q.tolist(),radius_m=e['radius_m']))
                counts['active']+=1
        return dict(physical_time_s=float(t),counts=dict(counts),active=active,completed_ids=finished,
            source_scene_sha256=self.sha256,pending_drawn_ids=[],display_interpolation_only=True)

    def partial_path(self,pid,elapsed):
        a=self.arrays[pid];position,_=self.position(pid,elapsed)
        k=np.searchsorted(a[:,0],elapsed,side='right')
        return np.vstack([a[:k,1:4],position]) if not np.array_equal(a[k-1,1:4],position) else a[:k,1:4]


def frame_schedule(scene,kind):
    """Physical timestamps explicit: sparse event windows or acquisition compression."""
    c=scene.catalog
    if kind=='02_accumulation':
        return [dict(time_s=float(t),label='ACQUISITION TIME COMPRESSION; completed tracks only',only_ids=None)
                for t in np.linspace(0,c['replay_end_time_s'],240)]
    if kind=='03_single_journeys':
        frames=[]
        for pid in sorted(c['representative_ids']):
            e=scene.entries[pid]
            for t in np.linspace(e['birth_time_s'],e['exit_time_s'],100):
                frames.append(dict(time_s=float(t),label=f'FULL JOURNEY MB {pid}; {e["exit_outlet"]}; SLOW MOTION; idle intervals CUT',only_ids=[pid]))
            for _ in range(12):
                frames.append(dict(time_s=float(e['exit_time_s']),label=f'OUTLET REACHED; MB {pid}; endpoint DISPLAY HOLD',only_ids=[pid]))
        return frames
    # Actual acquisition windows; rare arrivals are not multiplied into a cloud.
    reps=sorted(c['representative_ids'])
    frames=[]
    for pid in reps:
        e=scene.entries[pid]
        for t in np.linspace(e['birth_time_s'],e['exit_time_s'],80):
            frames.append(dict(time_s=float(t),label=f'EVENT WINDOW MB {pid}; SLOW MOTION; inter-window idle time CUT',only_ids=None))
    for t in np.repeat(c['replay_end_time_s'],20):
        frames.append(dict(time_s=float(t),label='CUT TO FINAL ACQUISITION ACCOUNTING',only_ids=None))
    return frames

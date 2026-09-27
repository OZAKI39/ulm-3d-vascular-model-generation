"""Read-only P7 replay. No solver, sampler, RNG, admission or geometry mutations.

SI is authoritative. Frames are right-continuous: admit -> active, exit -> deleted.
The ordered event ledger retains the zero-duration admitted/exited transitions.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import math
import numpy as np

REPO = Path(__file__).resolve().parents[3]
P7 = REPO / 'particle_3d/reports/particle7'
DEFAULT_OUTPUT = REPO / 'particle_3d/reports/particle8'
SCHEMA = 'PARTICLE8_REPLAY_V1'
ASSUMPTIONS = dict(H_D=0.45, C_MB=8.5e12, cumulative_flux_deterministic=True,
                   isotropic_orientation_v0=True, orientation_role='MODEL ASSUMPTION',
                   concentration_role='NOMINAL_ESTIMATE_NOT_MOUSE_SPECIFIC_MEASUREMENT',
                   H_D_role='FEED_VOLUME_FRACTION_NOT_INSTANTANEOUS_LUMEN_SNAPSHOT')
LIMITATIONS = [
    'REAL_RBC_CONTINUOUS_ADMISSION_AND_PASSAGE_NOT_ESTABLISHED',
    'SYNTHETIC_CONTROL_IS_NOT_PHYSIOLOGICAL_LUMEN_PROOF',
    'REAL_MIXED_SMOKE_HAS_STATIC_ADMITTED_OBSTACLES_NOT_MIXED_TRANSPORT',
    'SHORT_OPEN_SECTION_CENTER_PATH_IS_NOT_FULL_SHAPE_PASSAGE',
    'NONSPHERICAL_LUBRICATION_NOT_FROZEN', 'PRODUCTION_TIMESTEP_NOT_FROZEN',
    'PRODUCTION_NEIGHBOR_CUTOFF_AND_SKIN_NOT_FROZEN',
    'P65_HANDOFF_EVENT_TIME_CONVERGENCE_NOT_ESTABLISHED',
    'FULL_SUSPENSION_NOT_ESTABLISHED', 'NO_FULL_PK', 'NO_CFD', 'NO_PARTICLE_7_5',
    'SPHERE_NORMAL_NEARFIELD_ONLY_NO_FULL_MANY_BODY_MOBILITY',
    'NO_TANGENTIAL_OR_ROTATIONAL_LUBRICATION_COUPLING',
    'NO_GLYCOCALYX_ADHESION_ROUGHNESS_OR_SHELL_MECHANICS',
    'LAMMPS_STORAGE_NEIGHBORS_RESTART_ONLY_SINGLE_MPI_RANK',
    'ISOTROPIC_ORIENTATION_MODEL_ASSUMPTION_NOT_MEASURED',
    'FEED_HCT_NOT_FORCED_TUBE_HCT', 'MB_CONCENTRATION_NOMINAL_NOT_MEASURED',
    'DISPLAY_INTERPOLATION_IS_NOT_NEW_PHYSICS_OR_TIMESTEP_CONVERGENCE',
]


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def geometry(event):
    if event['species'] == 'MB':
        return dict(original_sample=deepcopy(event), radius_m=event['radius_m'],
                    diameter_m=event['diameter_um'] * 1e-6, volume_m3=event['volume_m3'],
                    orientation_wxyz=event['q'])
    g = event['geometry']
    return dict(original_sample=deepcopy(event), axes_m=[g[k] for k in ['a_m','b_m','c_m']],
                volume_m3=event['volume_m3'], aspect_ratio=g['c_m']/g['a_m'],
                orientation_wxyz=event['q'])


def from_population(raw, name, source_classification, horizon, *, restart_time=None):
    exits = {x['particle_id']: x for x in raw['exits']}
    active = {x['particle_id']: x for x in raw['active']}
    records = []
    for event in raw['events']:
        pid = event['particle_id']; birth = raw['births'].get(str(pid)); end = exits.get(pid)
        admitted = None if birth is None else birth['admitted_time_s']
        exit_time = None if end is None else end['exit_time']
        g = geometry(event)
        r = dict(particle_id=pid, species=event['species'], source_classification=source_classification,
                 geometry=g, scheduled_time_s=event['scheduled_time_s'], birth_time_s=admitted,
                 admit_time_s=admitted, exit_time_s=exit_time, delete_time_s=exit_time,
                 pending_duration_s=(horizon if admitted is None else admitted)-event['scheduled_time_s'],
                 pending_duration_censored=admitted is None, trajectory=[], admitted_shape=None,
                 outlet_role=None if end is None else end['outlet_role'])
        if birth is not None:
            p0 = birth['birth_position_m']; q = event['q']
            if source_classification == 'synthetic_control':
                shape = dict(mode='SPHERE_MB' if event['species']=='MB' else 'FREE_OBLATE',
                             radius_m=g.get('radius_m'), axes_m=g.get('axes_m'), q=q)
                last_time = horizon if end is None else exit_time
                last_position = active[pid]['position'] if end is None else end['trajectory_summary']['exit_position_m']
                r['trajectory'] = [dict(time_s=admitted, position_m=p0, q=q),
                                   dict(time_s=last_time, position_m=last_position, q=q)]
                r['trajectory_role'] = 'EXACT_COMMON_PLUG_FROM_P7_BIRTH_AND_EXIT_OR_FINAL_ENDPOINTS'
            else:
                # P7 mixed smoke deliberately holds admitted shapes static. Its stored
                # diagnostic velocity/omega must NOT be used as a propagation law.
                shape = deepcopy(active[pid]); shape['mode'] = birth['birth_shape_mode']
                r['trajectory'] = [dict(time_s=admitted,position_m=p0,q=q),
                                   dict(time_s=horizon,position_m=p0,q=q)]
                r['trajectory_role'] = 'STATIC_ADMISSION_OBSTACLE_NOT_TRANSPORT'
            r['admitted_shape'] = shape
        records.append(r)
    return make_scene(name, source_classification, [0.,horizon], records, restart_time)


def make_scene(name, source_classification, time_range, records, restart_time=None):
    ledger = []
    for r in records:
        for status,key in [('scheduled','scheduled_time_s'),('pending','scheduled_time_s'),
                           ('admitted','admit_time_s'),('active','admit_time_s'),
                           ('exited','exit_time_s'),('deleted','delete_time_s')]:
            if r[key] is not None:
                ledger.append(dict(time_s=r[key],particle_id=r['particle_id'],species=r['species'],status=status))
    order = {s:i for i,s in enumerate(['scheduled','pending','admitted','active','exited','deleted'])}
    ledger.sort(key=lambda e:(e['time_s'],e['particle_id'],order[e['status']]))
    scene = dict(schema=SCHEMA,name=name,source_classification=source_classification,
                 time_range_s=time_range,time_role='PHYSICAL_SIMULATION_CLOCK',units=dict(position='m',time='s',volume='m3'),
                 assumptions=ASSUMPTIONS,carried_forward_limitations=LIMITATIONS,
                 restart_time_s=restart_time,records=records,event_ledger=ledger,
                 interpolation='LINEAR_POSITION_SHORTEST_ARC_SLERP_DISPLAY_ONLY_NO_EXTRAPOLATION',
                 pending_position_policy='NULL_NO_INVENTED_UPSTREAM_SPATIAL_QUEUE',
                 birth_time_definition='ACTUAL_ADMISSION; SCHEDULED_TIME_STORED_SEPARATELY',
                 deletion_time_provenance='SAME_P7_LIFECYCLE_TRANSACTION_AS_EXIT')
    validate_scene(scene)
    return scene


def validate_scene(scene):
    if scene['schema'] != SCHEMA or scene['source_classification'] not in ['real_frozen','synthetic_control']:
        raise ValueError('Unknown replay schema or source classification')
    lo,hi = scene['time_range_s']
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi): raise ValueError('Invalid time range')
    ids = [r['particle_id'] for r in scene['records']]
    if len(ids)!=len(set(ids)) or any(not isinstance(i,int) or i<=0 for i in ids): raise ValueError('Duplicate/invalid ID')
    for r in scene['records']:
        if r['source_classification']!=scene['source_classification']: raise ValueError('Mixed source classes')
        if r['species'] not in ['MB','RBC']: raise ValueError('Unknown species')
        times = [r[k] for k in ['scheduled_time_s','admit_time_s','exit_time_s','delete_time_s'] if r[k] is not None]
        if not all(math.isfinite(t) and lo<=t<=hi for t in times) or times!=sorted(times): raise ValueError('Invalid event order/time')
        if (r['exit_time_s'] is not None and r['admit_time_s'] is None) or r['exit_time_s']!=r['delete_time_s']: raise ValueError('Invalid lifecycle')
        if r['geometry']['volume_m3']<=0: raise ValueError('Invalid geometry volume')
        if r['admit_time_s'] is None and r['trajectory']: raise ValueError('Pending particle cannot have a lumen trajectory')
        t = [s['time_s'] for s in r['trajectory']]
        if t and (len(t)<2 or np.any(np.diff(t)<=0) or t[0]!=r['admit_time_s'] or t[-1]!=(r['exit_time_s'] if r['exit_time_s'] is not None else hi)):
            raise ValueError('Incomplete trajectory interval')
        for s in r['trajectory']:
            if not np.isfinite(s['position_m']).all() or np.shape(s['position_m'])!=(3,): raise ValueError('Invalid position')
            if np.shape(s['q'])!=(4,) or not np.isclose(np.linalg.norm(s['q']),1,atol=1e-12,rtol=0): raise ValueError('Invalid quaternion')
    return True


def slerp(q0,q1,alpha):
    q0,q1 = np.asarray(q0),np.asarray(q1)
    if np.array_equal(q0,q1): return q0.copy()
    if alpha==0: return q0.copy()
    if alpha==1: return q1.copy()
    dot = float(q0@q1)
    if dot<0: q1=-q1; dot=-dot
    if dot>0.9995: q=(1-alpha)*q0+alpha*q1; return q/np.linalg.norm(q)
    theta=np.arccos(np.clip(dot,-1,1))
    return (np.sin((1-alpha)*theta)*q0+np.sin(alpha*theta)*q1)/np.sin(theta)


def interpolate(record,t):
    states=record['trajectory']
    if not states or t<states[0]['time_s'] or t>states[-1]['time_s']: raise ValueError('No extrapolation or pending trajectory')
    times=[s['time_s'] for s in states]
    i=min(int(np.searchsorted(times,t,side='right'))-1,len(states)-2)
    a,b=states[i:i+2]; f=(t-a['time_s'])/(b['time_s']-a['time_s'])
    # Preserve raw endpoints bit for bit; interpolate only display frames.
    p=a['position_m'] if f==0 or a['position_m']==b['position_m'] else b['position_m'] if f==1 else (np.array(a['position_m'])+f*(np.array(b['position_m'])-a['position_m'])).tolist()
    return dict(position_m=list(p),q=slerp(a['q'],b['q'],f).tolist())


def snapshot(scene,t):
    lo,hi=scene['time_range_s']
    if not math.isfinite(t) or not lo<=t<=hi: raise ValueError('Clock outside recorded scene')
    counts={sp:{k:0 for k in ['scheduled','admitted','pending','active','exited','deleted']} for sp in ['RBC','MB']}
    active=[]; pending=[]
    for r in scene['records']:
        if t<r['scheduled_time_s']: continue
        c=counts[r['species']]; c['scheduled']+=1
        if r['admit_time_s'] is None or t<r['admit_time_s']:
            c['pending']+=1
            pending.append(dict(particle_id=r['particle_id'],species=r['species'],status='pending',position_m=None,
                                pending_duration_s=t-r['scheduled_time_s']))
            continue
        c['admitted']+=1
        if r['exit_time_s'] is not None and t>=r['exit_time_s']: c['exited']+=1; c['deleted']+=1
        else:
            c['active']+=1
            active.append(dict(particle_id=r['particle_id'],species=r['species'],status='active',**interpolate(r,t)))
    for c in counts.values():
        assert c['scheduled']==c['pending']+c['admitted']
        assert c['admitted']==c['active']+c['exited'] and c['deleted']==c['exited']
    cut=scene['restart_time_s']
    return dict(time_s=t,source_classification=scene['source_classification'],counts=counts,
                active=active,pending=pending,restart_segment=('pre_restart' if t<cut else 'post_restart') if cut is not None else 'not_applicable')


def build_bundle(output=DEFAULT_OUTPUT):
    output=Path(output); folder=output/'data'; folder.mkdir(parents=True,exist_ok=True)
    used={}
    def load(name):
        p=P7/'data'/name; used[str(p.relative_to(REPO))]=digest(p); return read(p)
    upstream={}
    for stage in ['particle6_5','particle7']:
        tag=stage.replace('particle','PARTICLE')
        for ext in ['VALIDATION.json','REVIEW.md']:
            p=REPO/'particle_3d/reports'/stage/(tag+'_'+ext)
            upstream[str(p.relative_to(REPO))]=digest(p)
            p.read_text()  # Preserve and require both human and machine upstream evidence.
    validation=read(P7/'PARTICLE7_VALIDATION.json')
    for path,h in validation['source_sha256'].items():
        if digest(REPO/path)!=h: raise ValueError('Upstream scientific source changed: '+path)
    for path,h in validation['data_sha256'].items():
        if digest(P7/path)!=h: raise ValueError('Upstream P7 evidence changed: '+path)
    raw=load('10_continuous_final.json'); other=load('10_restarted_final.json')
    fields=['events','births','exits','pending','active','scheduler']
    parity={k:raw[k]==other[k] for k in fields}
    if not all(parity.values()): raise ValueError('Upstream restart mismatch')
    scenes={}
    scenes['synthetic']=from_population(raw,'synthetic_open_section','synthetic_control',raw['time_s'],restart_time=.125)
    scenes['restarted']=from_population(other,'synthetic_restarted','synthetic_control',other['time_s'],restart_time=.125)
    for k in ['synthetic','restarted']:
        scenes[k]['control_geometry']=dict(width_m=200e-6,length_m=1e-8,velocity_m_s=[0,0,2.5e-5],
            display_role='CENTER_PATH_SCHEMATIC_AXIAL_COORDINATE_STRETCHED_GLYPHS_NOT_TO_SCALE',
            finite_shapes_cross_open_caps=True,volume_definition='CENTER_ASSIGNED_WHOLE_PARTICLE')
    smoke=load('13_real_smoke.json')
    scenes['real_mixed']=from_population(smoke,'real_mixed_inlet_smoke','real_frozen',smoke['horizon_s'])
    mb=load('13_isolated_mb_transport.json'); event=mb['geometry']; admitted=event['scheduled_time_s']
    record=dict(particle_id=event['particle_id'],species='MB',source_classification='real_frozen',geometry=geometry(event),
                scheduled_time_s=admitted,birth_time_s=admitted,admit_time_s=admitted,exit_time_s=None,delete_time_s=None,
                pending_duration_s=0.,pending_duration_censored=False,outlet_role=None,
                admitted_shape=dict(mode='SPHERE_MB',radius_m=event['radius_m'],q=event['q']),
                trajectory=[dict(time_s=x['time_s'],position_m=x['particle']['position'],q=x['particle']['q']) for x in mb['states']],
                trajectory_role='P7_SAVED_P65_ACCEPTED_STATES_LINEAR_DISPLAY_INTERPOLATION')
    scenes['real_single_mb']=make_scene('real_single_mb_entry','real_frozen',[admitted-0.0001,mb['states'][-1]['time_s']],[record])
    queue=load('07_pending_queue.json'); original=queue['original']; birth=queue['admitted']
    pending_record=dict(particle_id=original['particle_id'],species='RBC',source_classification='synthetic_control',
        geometry=geometry(original),scheduled_time_s=original['scheduled_time_s'],birth_time_s=birth['admitted_time_s'],
        admit_time_s=birth['admitted_time_s'],exit_time_s=None,delete_time_s=None,
        pending_duration_s=birth['admitted_time_s']-original['scheduled_time_s'],pending_duration_censored=False,
        trajectory=[],admitted_shape=None,outlet_role=None)
    # Queue fixture ends at the instant of admission: preserve it as an event audit,
    # not as an invented continuation trajectory.
    write(folder/'pending_identity.json',dict(source_classification='synthetic_control',original=original,
        blocked=queue['blocked'],admitted=birth,times_s=queue['times_s'],queue_sizes=queue['queue_sizes'],
        identity_exact=all(original[k]==birth[k] and all(b[k]==original[k] for b in queue['blocked'])
                           for k in ['particle_id','geometry','q','provenance','scheduled_time_s']),
        replay_record=pending_record))
    cached=['01_inlet_mesh.json','01_flux_audit.json','02_sampling_statistics.json','13_real_geometry.json',
            '10_restart_parity.json','08_timestep_parity.json','14_tube_hct_statistics.json','11_independent_ledger_audit.json']
    for n in cached: write(folder/n,load(n))
    for n in ['02_actual_inlet_samples.csv','02_actual_triangle_counts.csv','03_04_scheduler.csv']:
        p=P7/'data'/n; used[str(p.relative_to(REPO))]=digest(p)
        (folder/n).write_bytes(p.read_bytes())
    for key,scene in scenes.items():
        scene['upstream_source_sha256']=used.copy(); write(folder/(key+'_scene.json'),scene)
    write(folder/'restart_exact_fields.json',dict(exact_match_fields=parity,checkpoint_time_s=.125,
        original_checkpoint_event_count=1180,final_event_count=len(raw['events']),
        full_rng_states_exact=raw['scheduler']['source']==other['scheduler']['source'],
        actual_destroy_and_binary_restore=load('10_restart_parity.json')['destroyed_actual_lammps_before_restart']))
    # Recompute final counters independently from time predicates, then compare P7.
    for key,counts in [('real_mixed',smoke['accounting']),('synthetic',raw['timeline'][-1])]:
        final=snapshot(scenes[key],scenes[key]['time_range_s'][1])
        for sp,c in final['counts'].items():
            for k in ['scheduled','admitted','pending','active','exited']:
                if c[k]!=counts[k+'_'+sp.lower()+'_count']: raise ValueError('Replay / P7 count mismatch')
    write(folder/'upstream_provenance.json',dict(upstream_reports_sha256=upstream,consumed_data_sha256=used,
        particle7_source_sha256=validation['source_sha256'],requested_mnt_data_files_present=False,
        resolution='REPOSITORY_COMMITTED_UPSTREAM_EQUIVALENTS',accepted_basis='USER_PARTICLE8_REQUEST_ACCEPTS_P65_AND_P7',
        p7_tested_commit=validation['git_commit'],assumptions=ASSUMPTIONS,
        source_files_verified=len(validation['source_sha256']),p7_data_files_verified=len(validation['data_sha256'])))
    return scenes

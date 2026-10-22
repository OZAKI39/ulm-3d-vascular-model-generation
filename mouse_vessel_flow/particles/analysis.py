"""Read-only reconstruction from original events, completion markers and raw samples."""
from collections import Counter
import json
import numpy as np
from ..paths import project_path, case_config
from ..provenance import sha256
from .formal_cohort_p9a5 import content_sha

def verify_tracks(name):
    cfg=case_config(name);spec=cfg.get('particles')
    if not spec:return dict(status='NOT_RUN',reason='No formal campaign attached to this baseline')
    cohort=json.loads(project_path(spec['cohort']).read_text());events=cohort['events']
    if len(events)!=spec['count']:raise ValueError('Cohort count changed')
    if [e['particle_id'] for e in events]!=list(range(1,len(events)+1)):raise ValueError('Particle source order changed')
    if any(b['source_event_id']<=a['source_event_id'] for a,b in zip(events,events[1:])):raise ValueError('Source events not ordered')
    outlets=Counter({o:0 for o in ('O1','O2','O3')});statuses=Counter();samples=0;hashes=0;maxdt=0.
    radius=[];safety=Counter();identities=[]
    for event in events:
        folder=project_path(spec['tracks'])/f"mb_{event['particle_id']:06d}"
        marker=json.loads((folder/'COMPLETE.json').read_text())
        assert marker['event_sha256']==content_sha(event)
        identity=marker['identity'];assert identity['flow_sha256']==cfg['flow']['sha256']
        if 'dt_s' in identity:assert identity['dt_s']==spec['dt_s']
        else:
            assert identity['formal_contract_sha256']==sha256(project_path(spec['contract']))
            contract=json.loads(project_path(spec['contract']).read_text())
            assert contract['dt_s']==spec['dt_s']
        if not identities:identities.append(identity)
        assert identity==identities[0]
        for relative,expected in marker['files'].items():
            path=(folder/relative).resolve()
            if not path.is_relative_to(folder.resolve()):raise ValueError('Completion path escapes track')
            assert sha256(path)==expected,(name,event['particle_id'],relative)
            hashes+=1
        meta=json.loads((folder/'trajectory.json').read_text());metrics=json.loads((folder/'metrics.json').read_text())
        with np.load(folder/'trajectory.npz',allow_pickle=False) as data:a=data['samples']
        assert a.ndim==2 and len(a)>0 and np.isfinite(a).all()
        assert np.all(np.diff(a[:,0])>0)
        if len(a)>1:maxdt=max(maxdt,float(np.diff(a[:,0]).max()))
        assert meta['integration_config']['dt_s']==spec['dt_s']
        assert np.array_equal(a[0,1:4],event['birth_center_m'])
        assert metrics['particle_id']==event['particle_id'] and metrics['source_event_id']==event['source_event_id']
        assert metrics['radius_m']==event['radius_m'] and metrics['diameter_um']==event['diameter_um']
        assert metrics['status']==marker['terminal_status']
        statuses[metrics['status']]+=1
        if metrics['outlet'] is not None:outlets[metrics['outlet']]+=1
        for key in ('penetration_count','handoff_violation_count','inlet_escape_count','nan_inf_count','unclassified_corruption_count'):
            safety[key]+=metrics[key]
        radius.append(event['radius_m']);samples+=len(a)
    assert dict(outlets)==spec['outlets'] and dict(statuses)==spec['statuses']
    assert maxdt<=spec['dt_s']+1e-12 and not any(safety.values())
    report=project_path(spec['report']);errors=report/'data/errors.json'
    if errors.exists():assert json.loads(errors.read_text())==[]
    return dict(status='PASS',count=len(events),outlets=dict(outlets),statuses=dict(statuses),
        failed_execution=0,raw_samples=samples,track_file_hashes_checked=hashes,
        dt_s=spec['dt_s'],maximum_accepted_dt_s=maxdt,source_event_last=events[-1]['source_event_id'],
        radius_range_m=[min(radius),max(radius)],safety_totals=dict(safety),original_identity=identities[0],
        interpretation='Read-only integrity/event checks; not a new timestep-convergence or CFD-accuracy proof')

def summarize(name):
    spec=case_config(name).get('particles')
    if not spec:return dict(case=name,status='NO_FORMAL_PARTICLE_CAMPAIGN')
    outlets=Counter({o:0 for o in ('O1','O2','O3')});statuses=Counter();count=0
    for f in sorted(project_path(spec['tracks']).glob('mb_*/metrics.json')):
        row=json.loads(f.read_text());count+=1;statuses[row['status']]+=1
        if row['outlet']:outlets[row['outlet']]+=1
    return dict(case=name,count=count,outlets=dict(outlets),statuses=dict(statuses),dt_s=spec['dt_s'])

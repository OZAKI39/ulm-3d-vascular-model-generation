"""P8.1 proposal/size/stop accounting by directly integrated point basins."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json
import numpy as np
from .particle82_provenance import atomic_json,sha256
from .particle82_integration import environment
from .injection_population import PopulationSource,FluxClock,LinearProfile,ConstantMBConcentrationV0,H_D,C_MB
from .particle7_cases import SONOVUE

OUTLETS=['OUTLET_01','OUTLET_02','OUTLET_03']


def csv_write(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    columns=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)


def prepare_p81_proposals(previous,output):
    previous=Path(previous);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    meta=[json.loads(p.read_text()) for p in sorted((previous/'trajectories').glob('*.json'))]
    positions=[];records=[];index={}
    for m in meta:
        indices=[]
        for candidate in m['admission_candidates']:
            i=len(positions);positions.append(candidate['position_m']);indices.append(i)
            records.append(dict(tracer_id=i,stable_id=m['particle_id'],draw=candidate['draw'],
                diameter_um=m['diameter_um'],status=candidate['status'],
                guard_exhausted=m['end_reason']=='INLET_ADMISSION_GUARD_EXHAUSTED_UNRESOLVED'))
        index[str(m['particle_id'])]=indices
    np.savez_compressed(output/'positions.npz',positions=np.array(positions))
    atomic_json(output/'PROPOSAL_INDEX.json',dict(records=records,indices_by_id=index,source='UNCHANGED_P81_ORIGINAL_PROPOSALS',
        original_metadata_sha256={p.name:sha256(p) for p in sorted((previous/'trajectories').glob('*.json'))}))
    return np.array(positions)


def basin_rows(point_output):
    folder=Path(point_output)
    if (folder/'shards').is_dir():folder=folder/'shards'
    files=sorted(folder.glob('tracer_*.json'))
    if not files:raise ValueError('No point-tracer shards in '+str(folder))
    rows=[]
    for p in files:rows.extend(json.loads(p.read_text())['rows'])
    rows.sort(key=lambda r:r['tracer_id'])
    if [r['tracer_id'] for r in rows]!=list(range(len(rows))):raise ValueError('Point basin IDs not complete/unique')
    return rows


def summarize_admission(previous,proposal_output,point_output):
    previous=Path(previous);out=Path(proposal_output)
    index=json.loads((out/'PROPOSAL_INDEX.json').read_text());points=basin_rows(point_output)
    if len(points)!=len(index['records']):raise ValueError('Proposal basin count mismatch')
    counter={b:Counter() for b in OUTLETS+['UNRESOLVED_POINT_PATH']}
    sizes={b:defaultdict(list) for b in counter};id_info={};csv_rows=[]
    for r,point in zip(index['records'],points):
        basin=point['outlet'] or 'UNRESOLVED_POINT_PATH';row=dict(r,basin=basin,point_end_reason=point['end_reason'])
        csv_rows.append(row);c=counter[basin];c['proposal_count']+=1
        if r['status']=='ACCEPTED':
            c['admitted_count']+=1;sizes[basin]['admitted_diameter_um'].append(r['diameter_um'])
        else:c['rejected_proposal_count']+=1
        if r['draw']==0:
            c['scheduled_first_proposal_count']+=1;sizes[basin]['scheduled_first_proposal_diameter_um'].append(r['diameter_um'])
            if r['guard_exhausted']:
                c['guard_exhausted_by_first_proposal_count']+=1;sizes[basin]['guard_exhausted_diameter_um'].append(r['diameter_um'])
    for pid,indices in index['indices_by_id'].items():
        rows=[csv_rows[i] for i in indices];accepted=[r for r in rows if r['status']=='ACCEPTED']
        id_info[pid]=dict(first_proposal_basin=rows[0]['basin'],admitted_basin=accepted[0]['basin'] if accepted else None,
            admission_attempts=len(rows),diameter_um=rows[0]['diameter_um'],guard_exhausted=rows[0]['guard_exhausted'])
    totals=dict(proposals=sum(c['proposal_count'] for c in counter.values()),
        scheduled=sum(c['scheduled_first_proposal_count'] for c in counter.values()),
        admitted=sum(c['admitted_count'] for c in counter.values()),
        guard_exhausted=sum(c['guard_exhausted_by_first_proposal_count'] for c in counter.values()))
    assert totals['scheduled']==2200 and totals['admitted']==1969 and totals['guard_exhausted']==231
    doc=dict(by_basin={b:dict(c) for b,c in counter.items()},diameters={b:dict(v) for b,v in sizes.items()},
        per_id=id_info,totals=totals,classification='DIRECT_ZERO_RADIUS_ADVECTION_OF_EACH_ORIGINAL_PROPOSAL_NOT_NEAREST_NEIGHBOR',
        accounting='Proposal counts include every retained retry. Scheduled/exhausted IDs use their first proposal basin; admitted IDs use their actually accepted position basin. Those distinct partitions are not interchangeable.',
        no_radius_resampling=True,unresolved_point_paths_preserved=True)
    atomic_json(out/'ADMISSION_BASIN_AUDIT.json',doc);csv_write(out/'proposal_by_basin.csv',csv_rows)
    return doc


def stop_audit(previous,admission,output):
    previous=Path(previous);out=Path(output);out.mkdir(parents=True,exist_ok=True);env=environment();rows=[]
    for p in sorted((previous/'trajectories').glob('*.json')):
        m=json.loads(p.read_text())
        if m['end_reason']!='INTEGRATION_SAFETY_STOP':continue
        a=np.load(p.with_suffix('.npz'))['samples'];last=a[-1];sample=env.field.sample(last[1:4])
        nominal=a[np.abs(a[:,0]/m['integration_config']['dt_s']-np.rint(a[:,0]/m['integration_config']['dt_s']))<1e-10]
        unchanged=0
        for k in range(len(nominal)-1,0,-1):
            if np.array_equal(nominal[k,1:4],nominal[k-1,1:4]):unchanged+=1
            else:break
        basin=admission['per_id'][str(m['particle_id'])]['admitted_basin'];recent=a[-64:]
        row=dict(stable_id=m['particle_id'],x_m=last[1],y_m=last[2],z_m=last[3],
            physical_time_s=m['birth_time_s']+last[0],elapsed_time_s=last[0],nearest_wall_gap_m=last[14],g_nf_m=last[15],h_lower_m=last[16],
            speed_m_s=np.linalg.norm(last[4:7]),nearfield_state_code=int(last[17]),tetra_id=sample.tetra_id,
            subdivision_depth=None,maximum_recorded_subdivision_depth=m['maximum_refinement_depth'],
            subdivision_depth_availability='P81_SAVED_MAXIMUM_ONLY_NOT_TERMINAL_TRIAL_DEPTH',provider_calls=m['provider_calls'],
            nominal_no_motion_count=unchanged,roundoff_motion_count=int(np.sum(np.linalg.norm(np.diff(recent[:,1:4],axis=0),axis=1)<=16*env.wall.roundoff_m)),
            last_accepted_displacement_m=np.linalg.norm(a[-1,1:4]-a[-2,1:4]),point_tracer_basin=basin,failure_detail=m['failure_detail'])
        rows.append({k:v.item() if isinstance(v,np.generic) else v for k,v in row.items()})
    assert len(rows)==720
    doc=dict(count=len(rows),rows=rows,by_basin=dict(Counter(r['point_tracer_basin'] for r in rows)),
        stop_is_computational_censoring_not_physiological_capture=True)
    atomic_json(out/'SAFETY_STOP_AUDIT.json',doc);csv_write(out/'safety_stop_records.csv',rows);return doc


def new_ledger(count=5000,seed=2026092182):
    env=environment();source=PopulationSource(SONOVUE,seed)
    clock=FluxClock(LinearProfile([0],[env.sampler.Q_m3_s]),ConstantMBConcentrationV0());events=[]
    for pid in range(1,count+1):
        event=source.next_mb();event['q']=source.orientation()[0].tolist()
        event.update(particle_id=pid,attempt_count=0,birth_time_s=clock.time_at(pid),scheduled_time_s=clock.time_at(pid),
                     position_seed=[seed,pid,82],orientation_role='ISOTROPIC_V0_MODEL_ASSUMPTION')
        events.append(event)
    return dict(events=events,seed=seed,source_rng_final=source.state(),Q_in_m3_s=env.sampler.Q_m3_s,
        C_MB_m3=C_MB,H_D_feed=H_D,acquisition_birth_window_s=clock.time_at(count),
        scheduler='P7_DETERMINISTIC_CUMULATIVE_FLUX',dataset_role='NATURAL_FLUX_WEIGHTED_DATASET',
        independent_superposition=True,no_diameter_resampling=True)

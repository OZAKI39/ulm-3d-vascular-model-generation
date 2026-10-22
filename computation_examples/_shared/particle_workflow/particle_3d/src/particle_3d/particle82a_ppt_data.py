"""User-revised scope: one 500-track presentation cohort, no population ranking."""
from pathlib import Path
from collections import Counter
import argparse,csv,json,shutil,socket,time
import numpy as np
from .particle82_provenance import atomic_json,sha256,require_remote
from .particle82a_pipeline import atomic_npz
from .particle82a_admission import context
from .particle82a_trajectories import initialize,run_tasks
from .particle8_replay import canonical_hash


def run(source,output,provenance,quota=500,workers=8):
    source=Path(source);output=Path(output);(output/'data').mkdir(parents=True,exist_ok=True)
    require_remote(provenance,provenance['hostname']);method='B';start=time.time()
    events=json.loads((source/'TRAJECTORY_COHORTS.json').read_text())['cohorts'][method][:quota]
    assert len(events)==quota
    atomic_json(output/'data/PPT_COHORT.json',dict(method=method,quota=quota,events=events,
        selection=f'FIRST_{quota}_ADMITTED_B_IDS_IN_COMMON_LEDGER_ORDER_NO_OUTLET_OR_OUTCOME_SELECTION',
        source_cohort_sha256=sha256(source/'TRAJECTORY_COHORTS.json'),
        role='PPT_DISPLAY_COHORT_NOT_LARGE_SAMPLE_STRATEGY_COMPARISON',
        size_role='FIXED_ANCHOR_POSITION_CONDITIONED_ACCEPTED_SONOVUE_SIZES',
        user_revision='STOP_15000_TRAJECTORY_STUDY; ONE_METHOD_ABOUT_500_TRACKS'))
    initialize(output,provenance);folder=output/'data/B/trajectories';folder.mkdir(parents=True,exist_ok=True)
    reused=[];tasks=[]
    for event in events:
        stem=f'mb_{event["particle_id"]:06d}'
        target=folder/(stem+'.json')
        if not target.exists():
            candidates=[source/'formal/B/trajectories'/(stem+'.json'),source/'benchmark/w8/B/trajectories'/(stem+'.json')]
            for old in candidates:
                if not old.exists():continue
                record=json.loads(old.read_text())
                if record['birth_metadata_sha256']!=canonical_hash(event):continue
                assert sha256(old.with_suffix('.npz'))==record['samples_sha256']
                for suffix in ['.json','.npz','.receipt.json']:shutil.copy2(old.with_suffix(suffix),folder/(stem+suffix))
                reused.append(dict(event_id=event['particle_id'],source=str(old),source_sha256=sha256(old)));break
        if target.exists():
            record=json.loads(target.read_text());receipt=json.loads(target.with_suffix('.receipt.json').read_text())
            assert record['birth_metadata_sha256']==canonical_hash(event)
            assert sha256(target.with_suffix('.npz'))==record['samples_sha256']==receipt['samples_sha256']
            assert sha256(target)==receipt['metadata_sha256']
        else:tasks.append((method,event,'data',.00025))
    performance={}
    if tasks:_,performance=run_tasks(tasks,workers,'PPT_500_TRACKS')
    c=context();records=[];checks=[];speeds=[];catalog=[]
    for event in events:
        file=folder/f'mb_{event["particle_id"]:06d}.json';record=json.loads(file.read_text());sample=file.with_suffix('.npz')
        receipt=json.loads(file.with_suffix('.receipt.json').read_text());a=np.load(sample)['samples']
        assert receipt['REMOTE_SERVER_COMPUTE'] and receipt['hostname']==socket.gethostname()
        assert sha256(sample)==record['samples_sha256'] and sha256(file)==receipt['metadata_sha256']
        assert len(a)>0 and np.isfinite(a).all() and np.all(np.diff(a[:,0])>0)
        assert np.array_equal(a[0,1:4],event['birth_center_m']) and event['birth_center_m']==event['anchor_m']
        if record['completed']:
            hit=c.env.classifier.first_event(a[-2,1:4],a[-1,1:4]);assert hit and hit.role==record['exit_outlet']
        speed=np.linalg.norm(a[:,4:7],axis=1)*1000;speeds.extend(speed.tolist())
        records.append(dict(event_id=event['particle_id'],diameter_um=record['diameter_um'],radius_m=record['radius_m'],
            birth_time_s=record['birth_time_s'],completed=record['completed'],outlet=record['exit_outlet'],
            end_reason=record['end_reason'],samples_path=str(sample.relative_to(output)),samples_sha256=sha256(sample),
            metadata_sha256=sha256(file),metadata_path=str(file.relative_to(output)),sample_count=len(a),
            last_age_s=float(a[-1,0]),residence_time_s=record['residence_time_s'],path_length_um=record['path_length_m']*1e6,
            computation_source_commit=receipt['git_commit'],compute_host=receipt['hostname']))
        catalog.append({k:records[-1][k] for k in ['event_id','diameter_um','completed','outlet','end_reason','sample_count','last_age_s','residence_time_s','path_length_um']})
    geometry={}
    for role,s in c.env.boundaries.items():geometry[role+'_points_m']=np.asarray(s.points);geometry[role+'_faces']=s.faces.reshape(-1,4)[:,1:]
    speed_nodes=np.linalg.norm(c.env.field.velocity_nodes_m_s,axis=1)
    geometry['branch_focus_m']=c.env.field.points[np.argmax(speed_nodes)] if hasattr(c.env.field,'points') else c.env.field.points_m[np.argmax(speed_nodes)]
    atomic_npz(output/'data/geometry.npz',**geometry)
    summary=dict(total_integrated=len(records),completed=sum(r['completed'] for r in records),
        outlets={b:sum(r['outlet']==b for r in records) for b in ['OUTLET_01','OUTLET_02','OUTLET_03']},
        end_reasons=dict(Counter(r['end_reason'] for r in records)),
        speed_max_mm_s=max(speeds),diameter_mean_um=float(np.mean([r['diameter_um'] for r in records])),
        diameter_range_um=[min(r['diameter_um'] for r in records),max(r['diameter_um'] for r in records)],
        original_sonovue_mean_um=c.distribution.target_mean_um(),
        frozen_inlet_mean_mm_s=c.env.sampler.Q_m3_s/c.geometry.area*1000)
    atomic_json(output/'data/SCENE.json',dict(records=records,summary=summary,
        cohort_sha256=sha256(output/'data/PPT_COHORT.json'),geometry_sha256=sha256(output/'data/geometry.npz'),
        frozen_input_sha256=provenance['frozen_input_sha256'],method='B',
        independent_tracks=True,position_interpolation='PIECEWISE_LINEAR_SAVED_ACCEPTED_SAMPLES',
        extrapolation=False,physical_radius_scale=1.,display_clock_role='INDEPENDENT_AGE_PHASE_OVERLAY_NOT_SIMULTANEOUS_PHYSIOLOGICAL_POPULATION'))
    with (output/'data/trajectory_catalog.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(catalog[0]));w.writeheader();w.writerows(catalog)
    atomic_json(output/'PPT_COMPUTE_VALIDATION.json',dict(all_pass=True,summary=summary,quota=quota,
        chosen_workers=workers,performance=performance,reused=reused,independently_computed_total=quota,
        REMOTE_SERVER_COMPUTE=True,hostname=socket.gethostname(),source_commit=provenance['source_git_commit'],
        all_sample_hashes_verified=True,all_actual_completed_outlet_segments_reclassified=True,
        frozen_input_sha256=provenance['frozen_input_sha256'],seconds=time.time()-start))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True)
    p.add_argument('--provenance',required=True);p.add_argument('--quota',type=int,default=500);p.add_argument('--workers',type=int,default=8)
    a=p.parse_args();run(a.source,a.output,json.loads(Path(a.provenance).read_text()),a.quota,a.workers)

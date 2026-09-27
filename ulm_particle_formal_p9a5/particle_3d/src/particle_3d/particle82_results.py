"""Export one provenance-bound natural scene; replay never samples new MBs."""
from pathlib import Path
from collections import Counter
import csv,gzip,json,os,socket,time
import numpy as np
from .particle8_replay import read,digest,canonical_hash
from .particle81_replay import Scene as PreviousScene,TITLE,FOOTERS,OUTLETS,INDEPENDENT
from .particle81_simulation import COLUMNS
from .particle82_provenance import atomic_json,require_remote
from .particle82_diagnostics import csv_write


def export_natural(root,host):
    from .particle82_tracers import init_environment,trace_positions
    root=Path(root);require_remote(host,host['hostname']);env=init_environment();started=time.time()
    data=root/'data';data.mkdir(exist_ok=True)
    run=read(root/'natural/catalog.json');births=read(root/'natural/birth_ledger.json')
    metas=[read(root/f'natural/trajectories/mb_{e["particle_id"]:06d}.json') for e in births['events']]
    positions=np.array([m.get('inlet_position_m',m['admission_candidates'][0]['position_m']) for m in metas])
    basins=trace_positions(positions,step_m=1e-7,error=1e-13)
    # Keep the exact point paths used to assign every natural ID's basin.
    offsets=np.cumsum([0]+[len(r['path']) for r in basins])
    np.savez_compressed(data/'natural_basin_point_paths.npz',positions=positions,offsets=offsets,
                        paths=np.concatenate([r['path'] for r in basins]))
    atomic_json(data/'NATURAL_BASIN_AUDIT.json',dict(hostname=socket.gethostname(),pid=os.getpid(),start_time=started,end_time=time.time(),
        git_commit=(Path(__file__).resolve().parents[3]/'SOURCE_COMMIT').read_text().strip(),
        frozen_input_sha256=host['frozen_input_sha256'],config=dict(step_m=1e-7,error=1e-13,horizon_m=2e-3),
        path_sha256=digest(data/'natural_basin_point_paths.npz'),role='DIAGNOSTIC_ONLY_ZERO_RADIUS_POSTPROCESS_OF_SAVED_NATURAL_ADMISSION_POSITIONS',
        rows=[dict(stable_id=m['particle_id'],outlet=r['outlet'],end_reason=r['end_reason']) for m,r in zip(metas,basins)]))
    entries=[];events=[];terminations=[];checks=[];geometry={}
    for role,s in env.boundaries.items():
        geometry[role+'_points_m']=np.asarray(s.points,dtype=float);geometry[role+'_faces']=s.faces.reshape(-1,4)[:,1:]
    np.savez_compressed(data/'full_frozen_geometry.npz',**geometry)
    atomic_json(data/'birth_ledger.json',births)
    with gzip.open(data/'trajectory_samples.csv.gz','wt',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['stable_id','absolute_time_s',*COLUMNS,'speed_m_s','active_until_terminal'])
        for m,point in zip(metas,basins):
            pid=m['particle_id'];meta_path=root/f'natural/trajectories/mb_{pid:06d}.json'
            sample_path=meta_path.with_suffix('.npz');receipt_path=root/f'natural/receipts/mb_{pid:06d}.json';receipt=read(receipt_path)
            if receipt['hostname']!=host['hostname'] or receipt['role']!='REMOTE_SERVER_HOST':raise ValueError('Local/wrong-host formal trajectory')
            if receipt['metadata_sha256']!=digest(meta_path) or m['samples_sha256']!=digest(sample_path):raise ValueError('Formal input changed')
            array=np.load(sample_path)['samples'];finite=bool(np.isfinite(array).all());monotone=bool(np.all(np.diff(array[:,0])>0))
            if not finite or not monotone:raise ValueError('Invalid formal samples')
            birth=births['events'][pid-1]
            if canonical_hash(birth)!=m['birth_metadata_sha256']:raise ValueError('Natural birth metadata changed')
            if len(array) and not np.array_equal(array[0,10:14],birth['q']):raise ValueError('Original quaternion changed')
            if m['completed']:
                hit=env.classifier.first_event(array[-2,1:4],array[-1,1:4])
                if hit is None or hit.role!=m['exit_outlet']:raise ValueError('Formal outlet reclassification failed')
            identity=float(np.max(np.abs(np.diff(array[:,1:4],axis=0)-np.diff(array[:,0])[:,None]*array[1:,4:7]))) if len(array)>1 else 0.
            checks.append(dict(stable_id=pid,finite=finite,monotone=monotone,terminal_reason=m['end_reason'],
                               birth_metadata_and_initial_orientation_verified=True,displacement_identity_max_error_m=identity,
                               official_outlet_verified=m['completed'],remote_host_verified=True))
            e={key:m.get(key) for key in ['particle_id','birth_time_s','radius_m','diameter_um','initial_q','inlet_face','inlet_barycentric',
                'inlet_position_m','completed','exit_outlet','exit_time_s','residence_time_s','path_length_m','end_reason','failure_detail',
                'sample_count','last_elapsed_time_s','last_physical_time_s','integration_config']}
            e.update(point_tracer_basin=point['outlet'] or 'UNRESOLVED_POINT_PATH',point_end_reason=point['end_reason'],
                point_seed_role='ACTUAL_ADMITTED_POSITION' if len(array) else 'FIRST_REJECTED_PROPOSAL_NO_ADMISSION',
                admission_attempts=len(m['admission_candidates']),admitted=bool(len(array)),
                samples_path=str(sample_path.relative_to(root)),metadata_path=str(meta_path.relative_to(root)),
                samples_sha256=digest(sample_path),metadata_sha256=digest(meta_path),
                receipt_path=str(receipt_path.relative_to(root)),receipt_sha256=digest(receipt_path),compute_host=receipt['hostname'],
                worker=receipt['worker_id'],config_sha256=receipt['config_sha256'],source_commit=receipt['git_commit'])
            entries.append(e)
            for i,row in enumerate(array):
                writer.writerow([pid,m['birth_time_s']+row[0],*row,float(np.linalg.norm(row[4:7])),int(i<len(array)-1)])
            events.append(dict(stable_id=pid,time_s=m['birth_time_s'],event='SCHEDULED',order=0))
            if len(array):
                events.append(dict(stable_id=pid,time_s=m['birth_time_s'],event='ADMITTED',order=1))
                events.append(dict(stable_id=pid,time_s=m['last_physical_time_s'],event=m['end_reason'],order=2))
                if m['completed']:events.append(dict(stable_id=pid,time_s=m['last_physical_time_s'],event='DELETED_AFTER_EXIT',order=3))
            else:events.append(dict(stable_id=pid,time_s=m['birth_time_s'],event=m['end_reason'],order=2))
            terminations.append(dict(stable_id=pid,birth_time_s=m['birth_time_s'],diameter_um=m['diameter_um'],radius_m=m['radius_m'],
                point_tracer_basin=e['point_tracer_basin'],admission_attempts=e['admission_attempts'],admitted=e['admitted'],
                termination=m['end_reason'],outlet=m['exit_outlet'],residence_time_s=m['residence_time_s'],path_length_m=m['path_length_m'],
                compute_host=receipt['hostname'],worker=receipt['worker_id'],config_sha256=receipt['config_sha256']))
    events.sort(key=lambda e:(e['time_s'],e['stable_id'],e['order']))
    completed=[e for e in entries if e['completed']];reps=[]
    for outlet in OUTLETS:
        cohort=sorted([e for e in completed if e['exit_outlet']==outlet],key=lambda e:(e['residence_time_s'],e['particle_id']))
        if cohort:reps.append(cohort[len(cohort)//2]['particle_id'])
    for e in sorted(completed,key=lambda e:e['path_length_m'],reverse=True):
        if len(reps)>=3:break
        if e['particle_id'] not in reps:reps.append(e['particle_id'])
    counts={o:sum(e['exit_outlet']==o for e in completed) for o in OUTLETS}
    scene=dict(schema='PARTICLE82_NATURAL_SAVED_SCENE_V1',title=TITLE,entries=entries,representative_ids=reps,
        units=dict(time='s',position='m',velocity='m/s',radius='m',diameter_um='um',gap='m',quaternion='dimensionless'),original_sample_columns=COLUMNS,
        scheduled=len(entries),admitted=sum(e['admitted'] for e in entries),completed=len(completed),outlet_counts=counts,
        end_reasons=dict(Counter(e['end_reason'] for e in entries)),acquisition_birth_window_s=births['acquisition_birth_window_s'],
        replay_end_time_s=max(e['last_physical_time_s'] or e['birth_time_s'] for e in entries),
        physical_samples=sum(e['sample_count'] for e in entries),dataset_role='NATURAL_FLUX_WEIGHTED_DATASET',
        independent_superposition=INDEPENDENT,scientific_labels=FOOTERS,formal_compute_host=host['hostname'],
        server_only_formal_dataset=True,source_commit=run['source_commit'],
        source_hashes=dict(birth_ledger=digest(data/'birth_ledger.json'),frozen_geometry=digest(data/'full_frozen_geometry.npz')),
        all_outlets_observed=all(counts.values()),full_network_visualization_target=all(n>=50 for n in counts.values()),
        display_only=dict(no_resampling=True,no_extrapolation=True,linear_position_interpolation=True,physical_coordinates_rotated=False))
    atomic_json(data/'particle8_2_natural_trajectory_catalog.json',scene)
    atomic_json(data/'trajectory_audit.json',dict(all_pass=True,rows=checks))
    atomic_json(data/'events.json',dict(events=events))
    csv_write(data/'events.csv',events);csv_write(data/'termination_records.csv',terminations)
    return scene


class Scene(PreviousScene):
    def __init__(self,root):
        self.root=Path(root);path=self.root/'data/particle8_2_natural_trajectory_catalog.json'
        self.catalog=read(path);self.sha256=digest(path)
        for key,name in [('birth_ledger','birth_ledger.json'),('frozen_geometry','full_frozen_geometry.npz')]:
            if digest(self.root/'data'/name)!=self.catalog['source_hashes'][key]:raise ValueError('Scene source changed')
        self.entries={e['particle_id']:e for e in self.catalog['entries']};self.arrays={}
        for pid,e in self.entries.items():
            for label in ['samples','metadata','receipt']:
                if digest(self.root/e[label+'_path'])!=e[label+'_sha256']:raise ValueError('Changed formal scene source')
            receipt=read(self.root/e['receipt_path'])
            if receipt['hostname']!=self.catalog['formal_compute_host']:raise ValueError('Formal mixed compute hosts')
            a=np.load(self.root/e['samples_path'])['samples'];a.flags.writeable=False;self.arrays[pid]=a
        self.completed=[pid for pid,e in self.entries.items() if e['completed']]

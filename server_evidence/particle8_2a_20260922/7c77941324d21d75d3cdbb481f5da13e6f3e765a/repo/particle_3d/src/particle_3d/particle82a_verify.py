"""Independent input, cohort, geometry and remote-execution evidence checks."""
from pathlib import Path
import argparse,gzip,json,socket,time
import numpy as np
from .particle82a_admission import context,common_event
from .particle82a_trajectories import iter_rows
from .particle82a_pipeline import CAUSES
from .particle82_provenance import atomic_json,sha256,require_remote
from .particle8_replay import canonical_hash


def verify(admission,trajectories,report,provenance):
    require_remote(provenance,provenance['hostname']);c=context();start=time.time()
    admission=Path(admission);trajectories=Path(trajectories);report=Path(report)
    counts=dict(events=0,original_checker_comparisons=0,A_B_trials=0,C_birth_geometry_checks=0,
                trajectory_records=0,trajectory_samples=0);max_birth_margin_error=0.;events_by_id={}
    for file in sorted((admission/'events').glob('events_*.json.gz')):
        stem=file.name.replace('.json.gz','');receipt=json.loads((file.parent/(stem+'.receipt.json')).read_text())
        assert receipt['hostname']==socket.gethostname() and receipt['REMOTE_SERVER_COMPUTE']
        for name,digest in receipt['files'].items():assert sha256(file.parent/name)==digest
        with gzip.open(file,'rt') as f:shard=json.load(f)
        raw=np.load(file.parent/(stem+'.npz'));attempts=raw['attempts']
        assert all(r['equal'] for r in shard['original_checker_parity'])
        counts['original_checker_comparisons']+=len(shard['original_checker_parity'])
        for row in shard['rows']:
            event=row['event'];pid=event['event_id'];counts['events']+=1
            original=common_event(pid)
            for key,value in original.items():assert event[key]==value,(pid,key)
            hashes={a['common_event_sha256'] for a in row['methods'].values()};assert len(hashes)==1
            assert next(iter(hashes))==__import__('hashlib').sha256(json.dumps(event,sort_keys=True).encode()).hexdigest()
            for mi,method in enumerate('AB'):
                result=row['methods'][method];trials=attempts[(attempts[:,0]==mi)&(attempts[:,1]==pid)]
                assert len(trials)==result['attempt_count'] and np.array_equal(trials[0,4:7],event['anchor_m'])
                assert trials[0,3]==event['first_diameter_um']
                if method=='A':assert np.all(trials[:,3]==event['first_diameter_um'])
                else:assert np.all(trials[:,4:7]==np.asarray(event['anchor_m']))
                assert bool(trials[-1,8]==0)==result['accepted']
                counts['A_B_trials']+=len(trials)
            result=row['methods']['C']
            assert result['anchor_m']==event['anchor_m'] and result['radius_m']==event['first_radius_m']
            assert not result['normal_fallback']
            if result['accepted']:
                margin=c.geometry.full_margin(result['birth_center_m'],result['radius_m'])
                assert margin>=0 and result['bracket_m'][1]-result['bracket_m'][0]<=result['minimality_tolerance_m']*1.000001
                max_birth_margin_error=max(max_birth_margin_error,abs(margin-result['full_margin_m']))
                counts['C_birth_geometry_checks']+=1
            events_by_id[pid]=row
    assert counts['events']>=100000
    assert len(events_by_id)==counts['events']
    assert sorted(events_by_id)==list(range(1,counts['events']+1))
    cohorts=json.loads((trajectories/'TRAJECTORY_COHORTS.json').read_text())['cohorts']
    for method in 'ABC':
        expected=[pid for pid,row in events_by_id.items() if row['methods'][method]['accepted']][:5000]
        assert [e['particle_id'] for e in cohorts[method]]==expected
        for event in cohorts[method]:
            pid=event['particle_id'];source=events_by_id[pid]['methods'][method]
            assert event['radius_m']==source['radius_m'] and event['birth_center_m']==source['birth_center_m']
            folder=trajectories/'formal'/method/'trajectories';file=folder/f'mb_{pid:06d}.json'
            record=json.loads(file.read_text());receipt=json.loads(file.with_suffix('.receipt.json').read_text())
            assert receipt['REMOTE_SERVER_COMPUTE'] and receipt['hostname']==socket.gethostname()
            assert receipt['metadata_sha256']==sha256(file) and record['samples_sha256']==sha256(file.with_suffix('.npz'))
            assert record['birth_metadata_sha256']==canonical_hash(event)
            samples=np.load(file.with_suffix('.npz'))['samples'];assert len(samples)==record['sample_count']>0
            assert np.array_equal(samples[0,1:4],source['birth_center_m'])
            assert np.isfinite(samples).all() and np.all(np.diff(samples[:,0])>=0)
            if record['completed']:
                hit=c.env.classifier.first_event(samples[-2,1:4],samples[-1,1:4])
                assert hit is not None and hit.role==record['exit_outlet']
            counts['trajectory_records']+=1;counts['trajectory_samples']+=len(samples)
    result=dict(all_pass=True,counts=counts,common_RNG_reproduction='EVERY_FORMAL_EVENT',
        first_A_B_candidate_parity=True,method_A_size_preserved=True,method_B_anchor_preserved=True,
        method_C_anchor_and_size_preserved=True,C_birth_margin_max_recheck_error_m=max_birth_margin_error,
        exact_cohort_selection_verified=True,all_actual_outlet_crossings_reclassified=True,
        all_saved_sample_hashes_verified=True,REMOTE_SERVER_COMPUTE=True,hostname=socket.gethostname(),
        source_commit=provenance['source_git_commit'],seconds=time.time()-start)
    atomic_json(report/'MACHINE_EVIDENCE_VERIFICATION.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--admission',required=True);p.add_argument('--trajectories',required=True)
    p.add_argument('--report',required=True);p.add_argument('--provenance',required=True);a=p.parse_args()
    print(json.dumps(verify(a.admission,a.trajectories,a.report,json.loads(Path(a.provenance).read_text())),indent=2))

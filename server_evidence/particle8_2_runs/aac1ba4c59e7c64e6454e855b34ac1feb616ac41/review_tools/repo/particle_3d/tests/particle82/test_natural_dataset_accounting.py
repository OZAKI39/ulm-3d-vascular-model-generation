import json,gzip,csv
from collections import Counter

def test_all_natural_ids_events_and_csv_states(saved_scene):
    s=saved_scene;c=s.catalog
    assert sorted(s.entries)==list(range(1,5001)) and c['scheduled']==5000
    assert c['admitted']==sum(bool(len(a)) for a in s.arrays.values())
    assert c['completed']==len(s.completed)==sum(c['outlet_counts'].values())
    assert c['dataset_role']=='NATURAL_FLUX_WEIGHTED_DATASET' and c['server_only_formal_dataset']
    assert c['all_outlets_observed']==all(c['outlet_counts'].values())
    final=s.snapshot(c['replay_end_time_s'])['counts'];assert final['active']==0
    assert final['scheduled']==final['completed']+final['stopped']+final['unresolved']
    rows=Counter()
    with gzip.open(s.root/'data/trajectory_samples.csv.gz','rt') as f:
        for row in csv.DictReader(f):rows[int(row['stable_id'])]+=1
    assert sum(rows.values())==c['physical_samples']
    assert all(rows[pid]==len(a) for pid,a in s.arrays.items())
    events=json.loads((s.root/'data/events.json').read_text())['events']
    assert events==sorted(events,key=lambda e:(e['time_s'],e['stable_id'],e['order']))

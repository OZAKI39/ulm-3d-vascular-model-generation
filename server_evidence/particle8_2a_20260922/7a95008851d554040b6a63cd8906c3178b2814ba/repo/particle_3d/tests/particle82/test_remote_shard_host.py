import json
from particle_3d.particle82_provenance import sha256

def test_all_formal_and_diagnostic_receipts(evidence_root,point_root):
    host=json.loads((evidence_root/'host_provenance.json').read_text())
    batches=[evidence_root/'natural',*(evidence_root/'scaling').glob('workers_[0-9]*'),*(evidence_root/'continuation').glob('guard_x*'),*(evidence_root/'timestep').glob('dt_div_*')]
    count=0
    for batch in batches:
        source=json.loads((batch/'catalog.json').read_text())['source_commit']
        for p in (batch/'receipts').glob('mb_*.json'):
            r=json.loads(p.read_text());meta=batch/'trajectories'/p.name
            assert r['hostname']==host['hostname']=='f7c62a262077'
            assert r['role']=='REMOTE_SERVER_HOST' and r['git_commit']==source
            assert r['metadata_sha256']==sha256(meta) and r['sample_sha256']==sha256(meta.with_suffix('.npz'))
            assert r['frozen_input_sha256']==host['frozen_input_sha256']
            assert r['start_time']<=r['end_time'];count+=1
    assert count>=5000+320+1440
    for folder in [point_root/'shards',evidence_root/'point_refinement_100000/shards',evidence_root/'admission/point_paths']:
        paths=list(folder.glob('tracer_*.json'));assert paths
        for p in paths:
            r=json.loads(p.read_text());assert r['hostname']==host['hostname']
            assert r['frozen_input_sha256']==host['frozen_input_sha256']
            assert r['path_sha256']==sha256(p.with_suffix('.npz'))

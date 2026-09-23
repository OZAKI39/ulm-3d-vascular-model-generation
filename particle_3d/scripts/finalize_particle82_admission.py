#!/usr/bin/env python3
"""Finalize retained, hash-verified proposal shards without reintegration."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_diagnostics import summarize_admission,stop_audit
from particle_3d.particle82_provenance import require_remote,sha256,atomic_json

p=argparse.ArgumentParser()
for n in ['previous','output-dir','host-provenance']:p.add_argument('--'+n,required=True)
a=p.parse_args();out=Path(a.output_dir);host=json.loads(Path(a.host_provenance).read_text())
require_remote(host,host['hostname'])
check=json.loads((out/'ORIGINAL_ADMISSION_RECHECK.json').read_text())
assert check['all_exact'] and check['checked_proposals']==149716
receipts=[]
for f in sorted((out/'point_paths').glob('tracer_*.json')):
    r=json.loads(f.read_text())
    assert r['hostname']==host['hostname'] and r['git_commit']==host['source_git_commit']
    assert r['path_sha256']==sha256(f.with_suffix('.npz'))
    receipts.append(dict(path=str(f),sha256=sha256(f)))
assert sum(json.loads(Path(r['path']).read_text())['count'] for r in receipts)==149716
admission=summarize_admission(a.previous,out,out/'point_paths')
stop=stop_audit(a.previous,admission,out/'stop_audit')
atomic_json(out/'AGGREGATION_RECEIPT.json',dict(source_commit=(Path(__file__).resolve().parents[2]/'SOURCE_COMMIT').read_text().strip(),
    integration_source_commit=host['source_git_commit'],retained_shards=receipts,reintegrated_proposals=0))
print(json.dumps(dict(admission=admission['by_basin'],stops=stop['by_basin'])),flush=True)

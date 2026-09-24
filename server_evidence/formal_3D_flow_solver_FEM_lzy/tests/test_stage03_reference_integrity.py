from stage03_helpers import *
from fem3d.audit import sha256

def test_all_frozen_history_and_original_shared_solver_sources_keep_their_bytes():
    baseline=read('reports/stage03/history_baseline.json')
    for name,record in baseline['files'].items():
        p=ROOT/name;assert p.is_file() and sha256(p)==record['sha256'],name
    for name,h in baseline['existing_source_sha256'].items():assert sha256(ROOT/name)==h,name
    c=config();assert sha256(Path(c['source']['configuration_path']))==c['source']['configuration_sha256']

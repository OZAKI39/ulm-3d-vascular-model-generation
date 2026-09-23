#!/usr/bin/env python3
"""Full read-only reference audit and frozen Stage 0/1 byte comparison."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import reference_snapshot,sha256,timestamp,write_json
baseline_path=ROOT/'reports/stage00/reference_baseline.json'
baseline=json.loads(baseline_path.read_text())
current=reference_snapshot(baseline['roots'])
write_json(ROOT/'reports/stage02/reference_final.json',current)
old,new=baseline['files'],current['files']
modified=sorted(k for k in old.keys()&new.keys() if old[k]!=new[k])
deleted=sorted(old.keys()-new.keys()); added=sorted(new.keys()-old.keys())
per_root={root:{name:sum(Path(p).is_relative_to(root) for p in paths) for name,paths in [('modified',modified),('deleted',deleted),('added',added)]} for root in baseline['roots']}
result={'status':'FAIL' if modified or deleted or added else 'PASS','timestamp':timestamp(),'entries':len(new),'per_root':per_root,'modified':modified,'deleted':deleted,'added':added,'baseline_sha256':sha256(baseline_path),'final_sha256':sha256(ROOT/'reports/stage02/reference_final.json'),'method':'Stage 0 full SHA256/size/mode/mtime inventory, including .git, caches and symlinks; atime excluded'}
write_json(ROOT/'reports/stage02/reference_integrity.json',result)
old=json.loads((ROOT/'reports/stage02/history_baseline.json').read_text())['files']; new={}
for category in ('reports','inputs','outputs','logs'):
    for stage in ('stage00','stage01'):
        for p in (ROOT/category/stage).rglob('*'):
            if p.is_file(): new[str(p.relative_to(ROOT))]={'sha256':sha256(p),'size':p.stat().st_size}
changes=sorted(k for k in old.keys()|new.keys() if old.get(k)!=new.get(k))
preservation={'status':'FAIL' if changes else 'PASS','timestamp':timestamp(),'entries':len(new),'changes':changes,'method':'Exact contents/size/file-set comparison to Stage 2 start; preserves the user formatting edit already present in reports/stage01/REPORT.md'}
write_json(ROOT/'reports/stage02/history_preservation.json',preservation)
print(json.dumps({'references':result,'frozen_history':preservation},indent=2))
sys.exit(int(result['status']!='PASS' or preservation['status']!='PASS'))

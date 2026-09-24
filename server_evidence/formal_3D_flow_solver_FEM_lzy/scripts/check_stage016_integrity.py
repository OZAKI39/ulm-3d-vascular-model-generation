#!/usr/bin/env python3
"""Audit all frozen stages and existing source without writing into them."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import reference_snapshot,sha256,timestamp,write_json
R=ROOT/'reports/stage01_6'
old_path=ROOT/'reports/stage00/reference_baseline.json';old=json.loads(old_path.read_text())
new=reference_snapshot(old['roots']);write_json(R/'reference_final.json',new)
a,b=old['files'],new['files'];modified=sorted(k for k in a.keys()&b.keys() if a[k]!=b[k]);deleted=sorted(a.keys()-b.keys());added=sorted(b.keys()-a.keys())
ref={'status':'FAIL' if modified or deleted or added else 'PASS','timestamp':timestamp(),'entries':len(b),'modified':modified,'deleted':deleted,'added':added,'per_root':{root:{name:sum(Path(p).is_relative_to(root) for p in paths) for name,paths in [('modified',modified),('deleted',deleted),('added',added)]} for root in old['roots']},'baseline_sha256':sha256(old_path),'final_sha256':sha256(R/'reference_final.json'),'method':'Full Stage 0 SHA256, size, mtime, mode, symlink and directory inventory; atime excluded'}
write_json(R/'reference_integrity.json',ref)
baseline=json.loads((R/'history_baseline.json').read_text());current={}
for category in ('reports','outputs','logs','inputs'):
    for stage in ('stage00','stage01','stage01_5','stage02'):
        for p in (ROOT/category/stage).rglob('*'):
            if p.is_file(): current[str(p.relative_to(ROOT))]={'sha256':sha256(p),'size':p.stat().st_size}
changes=sorted(k for k in baseline['files'].keys()|current.keys() if baseline['files'].get(k)!=current.get(k))
code_changes=[p for p,h in baseline['existing_source_sha256'].items() if not (ROOT/p).is_file() or sha256(ROOT/p)!=h]
history={'status':'FAIL' if changes or code_changes else 'PASS','timestamp':timestamp(),'file_count':len(current),'file_changes':changes,'existing_source_changes':code_changes,'method':'Frozen Stage 0/1/1.5/2 reports, outputs, logs and inputs: exact file set, SHA256 and size; all pre-existing fem3d source files including Stage 2 core hash matched; pre-existing user Stage 1 report formatting retained'}
write_json(R/'history_preservation.json',history)
print(json.dumps({'references':ref,'history':history},indent=2))
assert ref['status']==history['status']=='PASS'

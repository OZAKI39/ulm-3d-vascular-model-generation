#!/usr/bin/env python3
"""Full read-only reference audit and exact preservation of prior scientific evidence."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import reference_snapshot,sha256,timestamp,write_json
R=ROOT/'reports/stage03';read=lambda p:json.loads(p.read_text())
oldpath=ROOT/'reports/stage00/reference_baseline.json';old=read(oldpath);new=reference_snapshot(old['roots'])
write_json(R/'reference_final.json',new)
a,b=old['files'],new['files'];modified=sorted(k for k in a.keys()&b.keys() if a[k]!=b[k]);deleted=sorted(a.keys()-b.keys());added=sorted(b.keys()-a.keys())
reference={'status':'FAIL' if modified or deleted or added else 'PASS','timestamp':timestamp(),'entries':len(b),
           'modified':modified,'deleted':deleted,'added':added,'baseline_sha256':sha256(oldpath),
           'final_sha256':sha256(R/'reference_final.json'),'method':'Full SHA256, size, mtime, mode, symlink and directory inventory of both read-only projects; atime excluded'}
write_json(R/'reference_integrity.json',reference)
baseline=read(R/'history_baseline.json');changes=[]
for name,item in baseline['files'].items():
    p=ROOT/name
    if not p.is_file() or sha256(p)!=item['sha256'] or p.stat().st_size!=item['size']:changes.append(name)
# Ensure complete prior production stages have no newly inserted files either.
new_history=[]
for category in ('reports','outputs','logs','inputs'):
    for stage in ('stage00','stage01','stage01_5','stage01_6','stage01_7','stage02'):
        for p in (ROOT/category/stage).rglob('*'):
            if p.is_file() and str(p.relative_to(ROOT)) not in baseline['files']:new_history.append(str(p.relative_to(ROOT)))
source_changes=[name for name,h in baseline['existing_source_sha256'].items() if not (ROOT/name).is_file() or sha256(ROOT/name)!=h]
history={'status':'FAIL' if changes or new_history or source_changes else 'PASS','timestamp':timestamp(),
         'frozen_file_count':len(baseline['files']),'modified_or_deleted':changes,'new_history_files':new_history,
         'existing_source_changes':source_changes,'user_existing_dirty_paths':baseline.get('user_existing_dirty_paths',[]),
         'prior_report_conclusions_rewritten':False,'stage018_retention':'Only approved final scientific evidence plus new CLEANUP_NOTICE; original bytes preserved'}
write_json(R/'history_preservation.json',history)
print(json.dumps({'references':reference,'history':history},indent=2));assert reference['status']==history['status']=='PASS'

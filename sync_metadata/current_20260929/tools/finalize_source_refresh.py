"""Supplement the full hash check for files edited concurrently during copy."""
from pathlib import Path
import ast,csv,datetime,hashlib,importlib.util,json,shutil
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260929'
spec=importlib.util.spec_from_file_location('collector',M/'tools/collect_snapshot.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
v=json.loads((M/'validation/snapshot_validation.json').read_text())
assert not v['failures'] and not v['credential_pattern_findings'] and not v['oversized_files']
shutil.copy2(M/'validation/snapshot_validation.json',M/'validation/full_check_before_final_refresh.json')
c.ROWS={r['destination']:r for r in csv.DictReader((M/'source_files.csv').open())}
trans={r['path']:r for r in json.loads((M/'documentation_transforms.json').read_text())}
refresh=[];checked=0
for rel,row in list(c.ROWS.items()):
    source=Path(row['source']);target=D/rel
    assert source.exists(),str(source)
    actual=source.stat();expected=target.stat()
    if row['source'] not in v['changed_sources'] and rel not in trans and actual.st_size==expected.st_size and actual.st_mtime_ns==expected.st_mtime_ns:continue
    checked+=1;h=c.sha(source)
    if h==row['sha256']:continue
    assert rel.startswith('vascular_printing/'), 'Scientific source changed after full validation: '+rel
    c.copy(source,Path(rel));assert c.sha(target)==c.ROWS[rel]['sha256']==c.sha(source)
    if target.suffix=='.py':ast.parse(target.read_text())
    refresh.append(dict(path=rel,before_sha256=row['sha256'],after_sha256=c.ROWS[rel]['sha256']))
for row in c.ROWS.values():
    if row['source'] in v['changed_sources']:assert c.sha(Path(row['source']))==row['sha256']
with (M/'source_files.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['source','destination','bytes','sha256']);w.writeheader();w.writerows(c.ROWS.values())
v.update(PASS=True,changed_sources=[],final_refresh=refresh,followup_changed_timestamp_hash_checks=checked,
    full_check_record='full_check_before_final_refresh.json',snapshot_cutoff_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    validation_scope='Full SHA256 check followed by incremental SHA256 refresh of concurrently edited files; other development was not locked')
(M/'validation/snapshot_validation.json').write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')
summary=json.loads((M/'collection_summary.json').read_text());summary['logical_bytes']=sum(int(r['bytes']) for r in c.ROWS.values())
(M/'collection_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(dict(PASS=True,refresh=refresh,followup_checks=checked),ensure_ascii=False,indent=2))

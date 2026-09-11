"""One-time pre-upload addition of eight source/evidence files; existing frozen files are unchanged."""
from transfer import *

def main():
    d=directory('20260910T212439Z')
    assert read(d/'transfer_status.json')['status']=='PREPARED'
    assert not (d/'metadata/remote_before.json').exists(),'Cannot amend a started remote transfer'
    manifest=validate_frozen(d);old_identity=read(d/'.transfer_identity.json')
    items,excluded,special=inventory();old={x['path'] for x in manifest['files']+manifest['symlinks']}
    assert not (old-set(items)) and not special
    added=sorted(set(items)-old);assert len(added)==8,'Inspect unexpected selection difference'
    audit=scan_sensitive({p:items[p]['source_path'] for p in added})
    assert audit['status']=='PASS'
    evidence=[];consistency=read(d/'metadata/source_consistency.json')
    for key in added:
        p=Path(items[key]['source_path']);q=d/'source'/key;a=identity(p);h=sha(p)
        q.parent.mkdir(parents=True,exist_ok=True)
        record=command_log(RSYNC+[str(p),str(q.parent)+'/'],BASE/'runtime_logs'/(d.name+'-supplement-'+hashlib.sha256(key.encode()).hexdigest()[:12]+'.log'))
        assert record['exit_code']==0 and identity(p)==a and sha(p)==sha(q)==h==audit['file_sha256'][key]
        assert q.stat().st_mode&0o111==a['mode']&0o111
        entry=dict(path=key,source_path=str(p),size=q.stat().st_size,sha256=h,execute_bits=q.stat().st_mode&0o111,
            original_copy=True,explicit_supplemental=True,redacted=False,link_rewrite=False)
        manifest['files'].append(entry);evidence.append(entry)
        consistency['before'][key]=consistency['after'][key]=dict(identity=a,sha256=h,target=None)
        items[key]['explicit_supplemental']=True
        items[key]['supplement_reason']='Pre-upload correction: build name patterns apply only to directories, preserve build scripts/records'
    amendment=d/'metadata/pre_upload_amendment';amendment.mkdir()
    shutil.copy2(d/'transfer_manifest.json',amendment/'original_transfer_manifest.json')
    shutil.copy2(d/'SHA256SUMS',amendment/'original_SHA256SUMS')
    write_json(amendment/'record.json',dict(recorded_at=now(),reason='Correct directory-only build exclusions before first upload; add eight omitted scripts/records without regenerating or modifying existing source copies',
        original_identity=old_identity,added_files=evidence,existing_frozen_source_files_modified=False,remote_transfer_started=False))
    manifest['files'].sort(key=lambda x:x['path']);manifest['regular_files']=len(manifest['files']);manifest['regular_bytes']=sum(x['size'] for x in manifest['files'])
    manifest['pre_upload_supplement_record']='metadata/pre_upload_amendment/record.json'
    write_json(d/'transfer_manifest.json',manifest)
    (d/'SHA256SUMS').write_text(''.join(x['sha256']+'  '+x['path']+'\n' for x in manifest['files']))
    old_identity.update(manifest_sha256=sha(d/'transfer_manifest.json'),sums_sha256=sha(d/'SHA256SUMS'))
    write_json(d/'.transfer_identity.json',old_identity)
    consistency['pre_upload_additions']=added;write_json(d/'metadata/source_consistency.json',consistency)
    sensitive=read(d/'metadata/sensitive_information_check.json');sensitive['file_sha256'].update(audit['file_sha256'])
    sensitive['scanned_files_and_members']+=audit['scanned_files_and_members'];write_json(d/'metadata/sensitive_information_check.json',sensitive)
    write_json(d/'metadata/native_source_locations.json',native_sources(items))
    write_json(d/'metadata/excluded_files.json',dict(exclusions=excluded,pattern_exclusions_are_subtree_records=True,source_files_modified=False))
    write_json(d/'metadata/supplemented_files.json',[dict(path=k,**v) for k,v in items.items() if v['explicit_supplemental']])
    for p in (BASE/'tools').iterdir():
        if p.is_file():shutil.copy2(p,d/'metadata/transfer_tools'/p.name)
    rules=(d/'rsync-exclude.txt').read_text()
    for pattern in BUILD_DIR_PATTERNS:rules=rules.replace('\n'+pattern+'\n','\n'+pattern+'/\n')
    (d/'rsync-exclude.txt').write_text(rules)
    with (d/'metadata/transfer_log.txt').open('a') as f:f.write(now()+' Added eight build scripts/records before first upload; no original source or existing frozen source file modified.\n')
    state(d,'PREPARED',source_regular_files=manifest['regular_files'],source_regular_bytes=manifest['regular_bytes'])
    validate_frozen(d)
    say(json.dumps(dict(status='PREPARED',added=[x['path'] for x in evidence],statistics=packet_stats(d)),ensure_ascii=False))

if __name__=='__main__':main()

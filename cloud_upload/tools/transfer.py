#!/usr/bin/env python3
"""Frozen local source snapshots and bounded rsync/SSH delivery. Standard library only."""
import argparse
from datetime import datetime,timezone
import getpass
import platform
import shlex
import shutil
import sys
import time
import uuid
from scope import *
from evidence import native_sources,inputs

SSH=['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=15','-o','StrictHostKeyChecking=yes',
     '-o','ClearAllForwardings=yes','-o','RemoteCommand=none']
RSYNC=['rsync','-a','--no-owner','--no-group','--no-devices','--no-specials']
VERIFIED={'TRANSFER_VERIFIED','TRANSFER_VERIFIED_INPUTS_PARTIAL'}

def now():return datetime.now(timezone.utc).isoformat()
def say(message):print(message,flush=True)
def checked_id(value):
    if not re.fullmatch(r'\d{8}T\d{6}Z(?:-[a-f0-9]{8})?',value):raise ValueError('Invalid TRANSFER_ID')
    return value
def directory(value):
    p=BASE/checked_id(value)
    if not p.is_dir() or p.is_symlink():raise RuntimeError('Transfer directory missing or unsafe')
    return p
def read(p):return json.loads(Path(p).read_text())
def state(d,status,**kwargs):
    previous=read(d/'transfer_status.json') if (d/'transfer_status.json').exists() else {}
    previous.update(transfer_id=d.name,upload_dir=str(d),remote_dir='/workspace/bloodflow/uploads/'+d.name,status=status,updated_at=now(),**kwargs)
    write_json(d/'transfer_status.json',previous)

def local_check():
    r=subprocess.run(['ssh','-G','vast-mirheo'],capture_output=True,text=True,check=True)
    config={}
    for line in r.stdout.splitlines():
        parts=line.split(None,1)
        if len(parts)==2 and parts[0] in ('hostname','port','user','identityfile'):config.setdefault(parts[0],[]).append(parts[1])
    assert config['hostname']==['204.111.105.196'] and config['port']==['50271'] and config['user']==['root'],'SSH target mismatch'
    assert config['identityfile']==['~/.ssh/id_ed25519_vast'] or config['identityfile']==['/home/lzy/.ssh/id_ed25519_vast'],'SSH identity path mismatch'
    info=dict(user=getpass.getuser(),home=os.environ.get('HOME'),kernel=platform.release(),cwd=os.getcwd(),ssh_effective_config=config,
        host_key_validation='StrictHostKeyChecking=yes; existing known_hosts unchanged',port_forwarding=False,
        tools={x:shutil.which(x) for x in ['ssh','rsync','python3']})
    assert info['user']=='lzy' and info['home']=='/home/lzy' and 'WSL' in info['kernel'].upper(),'Must run in local lzy WSL'
    assert all(info['tools'].values()),'Required existing tool missing'
    return info

def command_log(argv,path,input_data=None):
    started=now();t=time.monotonic()
    with path.open('wb') as log:
        result=subprocess.run(argv,input=input_data,stdout=log,stderr=subprocess.STDOUT)
    record=dict(argv=argv,started_at=started,finished_at=now(),elapsed_s=time.monotonic()-t,exit_code=result.returncode,log_path=str(path))
    write_json(path.with_suffix(path.suffix+'.json'),record)
    return record

def remote(mode,d=None,**kwargs):
    request=dict(mode=mode,**kwargs)
    if d:
        request.update(transfer_id=d.name,remote_dir='/workspace/bloodflow/uploads/'+d.name,identity=read(d/'.transfer_identity.json'))
    # Request JSON travels as Python string data on stdin, never as a shell argument.
    payload=('import sys\nsys.argv = ["remote_check", '+repr(json.dumps(request))+']\n'+(BASE/'tools/remote_check.py').read_text()).encode()
    argv=SSH+['vast-mirheo',shlex.join(['python3','-B','-'])]
    r=subprocess.run(argv,input=payload,capture_output=True)
    lines=r.stdout.decode(errors='replace').splitlines()
    parsed=None
    for line in reversed(lines):
        try:parsed=json.loads(line);break
        except json.JSONDecodeError:pass
    if r.returncode or not parsed or parsed.get('status')!='PASS':
        # SSH diagnostics may contain hostnames/key paths, but never private key contents.
        failure=dict(mode=mode,ssh_exit_code=r.returncode,result=parsed,ssh_diagnostics=r.stderr.decode(errors='replace')[:2000])
        write_json(BASE/'runtime_logs'/((d.name if d else 'preflight')+'-remote-failure.json'),failure)
        raise RuntimeError('SSH_OR_REMOTE_CHECK_FAILED: '+json.dumps(failure,ensure_ascii=False))
    return parsed

def packet(d):
    paths=[]
    for root,dirs,files in os.walk(d,followlinks=False):
        for name in dirs+files:
            p=Path(root)/name
            if p.is_file() or p.is_symlink():paths.append(p.relative_to(d).as_posix())
    return sorted(paths)

def packet_stats(d):
    rows=[dict(path=relative,size=(d/relative).lstat().st_size) for relative in packet(d) if not (d/relative).is_symlink()]
    return dict(regular_files=len(rows),regular_bytes=sum(x['size'] for x in rows),symlinks=sum((d/p).is_symlink() for p in packet(d)),largest20=sorted(rows,key=lambda x:x['size'],reverse=True)[:20])

def freeze_files(d,items):
    before={};rows=[];links=[]
    for key,item in sorted(items.items()):
        if any(c in key for c in ['\n','\r','\\']):raise RuntimeError('Unsupported checksum filename: '+repr(key))
        p=Path(item['source_path']);ident=identity(p)
        before[key]=dict(identity=ident,sha256=sha(p) if p.is_file() and not p.is_symlink() else None,
            target=os.readlink(p) if p.is_symlink() else None)
        if identity(p)!=ident:raise RuntimeError('Source changed while hashing: '+key)
    for project in ROOTS:
        filelist=BASE/'runtime_logs'/(d.name+'-'+project+'-files.nul')
        filelist.write_bytes(b''.join((k[len(project)+1:].encode()+b'\0') for k in sorted(items) if k.startswith(project+'/')))
        target=d/'source'/project;target.mkdir(parents=True,exist_ok=True)
        log=command_log(RSYNC+['-r','--from0','--files-from='+str(filelist),str(ROOTS[project])+'/',str(target)+'/'],BASE/'runtime_logs'/(d.name+'-copy-'+project+'.log'))
        if log['exit_code']:raise RuntimeError('Local rsync copy failed: '+project)
    changed=[];after={}
    for key,item in sorted(items.items()):
        p=Path(item['source_path']);q=d/'source'/key;a=identity(p)
        h=sha(p) if p.is_file() and not p.is_symlink() else None
        after[key]=dict(identity=a,sha256=h,target=os.readlink(p) if p.is_symlink() else None)
        if after[key]!=before[key] or identity(p)!=a:changed.append(key)
    if changed:
        write_json(d/'metadata/source_changes.json',dict(paths=changed,retries=0))
        # A bounded retry refreshes only the affected files, with explicit evidence.
        for key in changed:
            p=Path(items[key]['source_path']);q=d/'source'/key;a=identity(p)
            h=sha(p) if not p.is_symlink() else None
            retry=command_log(RSYNC+['--checksum',str(p),str(q.parent)+'/'],BASE/'runtime_logs'/(d.name+'-retry-'+hashlib.sha256(key.encode()).hexdigest()[:12]+'.log'))
            if retry['exit_code'] or identity(p)!=a or (h is not None and sha(p)!=h):raise RuntimeError('Source still changing; stopped affected copy: '+key)
            after[key]=dict(identity=a,sha256=h,target=os.readlink(p) if p.is_symlink() else None)
        write_json(d/'metadata/source_changes.json',dict(paths=changed,retries=1,status='AFFECTED_FILES_RECOPIED_AND_STABLE'))
    for key,item in sorted(items.items()):
        p=Path(item['source_path']);q=d/'source'/key
        base=dict(path=key,source_path=str(p),original_copy=True,explicit_supplemental=item['explicit_supplemental'],redacted=False,link_rewrite=False)
        if p.is_symlink():
            if not q.is_symlink() or os.readlink(q)!=after[key]['target']:raise RuntimeError('Local link copy mismatch: '+key)
            links.append(dict(base,type='relative_internal_symlink',target=os.readlink(q)))
        else:
            s=q.stat()
            if sha(q)!=after[key]['sha256'] or s.st_size!=after[key]['identity']['size'] or s.st_mode&0o111!=after[key]['identity']['mode']&0o111:raise RuntimeError('Local copy hash/mode mismatch: '+key)
            rows.append(dict(base,size=s.st_size,sha256=sha(q),execute_bits=s.st_mode&0o111))
    current,_,_=inventory()
    if set(current)!=set(items):raise RuntimeError('Source selection changed during snapshot; pause preparation')
    write_json(d/'metadata/source_consistency.json',dict(status='PASS',before=before,after=after,recopied_files=changed,
        original_files_unchanged_during_copy=not changed,source_and_staging_sha256_and_execute_bits_equal=True,
        note='Filesystem snapshot only; unsaved editor buffers are not included. Concurrent process activity checked separately.'))
    return rows,links

def prepare():
    local=local_check()
    for old in sorted(BASE.iterdir()):
        if old.is_dir() and re.fullmatch(r'\d{8}T\d{6}Z(?:-[a-f0-9]{8})?',old.name):
            oldstatus=read(old/'transfer_status.json').get('status') if (old/'transfer_status.json').exists() else 'UNKNOWN'
            if oldstatus not in VERIFIED:raise RuntimeError('Unfinished transfer exists; inspect/resume same ID '+old.name+' status='+oldstatus)
    items,excluded,special=inventory();process=process_audit()
    if process['open_source_write_descriptors']:raise RuntimeError('Open write handles in source; no process was stopped')
    if special:raise RuntimeError('Unsupported special files in selection')
    linkcheck=links_audit(items)
    if linkcheck['status']!='PASS':raise RuntimeError('Unresolved selected symlink; affected copy stopped')
    say('Auditing selected files and archive members before freezing...')
    sensitive=scan_sensitive({k:v['source_path'] for k,v in items.items()})
    write_json(BASE/'runtime_logs/prepared_sensitive_audit.json',sensitive)
    if sensitive['status']!='PASS':raise RuntimeError('Sensitive information audit requires review; see runtime_logs/prepared_sensitive_audit.json')
    tid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    if (BASE/tid).exists():tid+='-'+uuid.uuid4().hex[:8]
    d=BASE/tid;d.mkdir();(d/'metadata').mkdir();(d/'source').mkdir()
    remote_dir='/workspace/bloodflow/uploads/'+tid
    (BASE/'current_transfer.env').write_text('\n'.join(k+'='+shlex.quote(v) for k,v in dict(UPLOAD_DIR=str(d),REMOTE_DIR=remote_dir,TRANSFER_ID=tid).items())+'\n')
    state(d,'TRANSFER_INCOMPLETE',phase='LOCAL_PREPARATION')
    try:
        native=native_sources(items);dependency,mapping,critical=inputs(items,remote_dir)
        write_json(d/'metadata/local_environment.json',local)
        write_json(d/'metadata/process_check_before.json',process)
        write_json(d/'metadata/project_git_status.json',[git_info(p) for p in ROOTS.values()])
        write_json(d/'metadata/native_source_locations.json',native)
        write_json(d/'metadata/input_dependency_check.json',dependency)
        write_json(d/'metadata/path_mapping.json',mapping)
        write_json(d/'metadata/excluded_files.json',dict(exclusions=excluded,pattern_exclusions_are_subtree_records=True,source_files_modified=False))
        write_json(d/'metadata/supplemented_files.json',[dict(path=k,**v) for k,v in items.items() if v['explicit_supplemental']])
        write_json(d/'metadata/sensitive_information_check.json',sensitive)
        write_json(d/'metadata/symlink_check.json',linkcheck)
        say('Freezing '+str(len(items))+' selected regular files and symlinks into '+str(d))
        files,links=freeze_files(d,items)
        if any(sensitive['file_sha256'].get(x['path'])!=x['sha256'] for x in files):raise RuntimeError('A file changed after its sensitive-information scan; stop before upload')
        # Recheck native provenance and selected-file identity after the copy.
        if native_sources(items)!=native:raise RuntimeError('Native source/provenance changed during preparation')
        write_json(d/'metadata/process_check_after.json',process_audit())
        write_json(d/'metadata/symlink_check.json',links_audit(items,d/'source'))
        manifest=dict(schema_version=1,transfer_id=tid,created_at=now(),source_roots={k:str(v) for k,v in ROOTS.items()},
            files=files,symlinks=links,regular_files=len(files),regular_bytes=sum(x['size'] for x in files),
            critical_files=critical,critical_directories=['mirheo_starter/py_scripts','mirheo_starter/scripts','mirheo_starter/test_code',
                'mirheo_starter/data/single_rbc_repair/'+CAMPAIGN+'/native/source','hemocell_starter/cases','hemocell_starter/vendor/HemoCell/palabos'])
        write_json(d/'transfer_manifest.json',manifest)
        (d/'SHA256SUMS').write_text(''.join(x['sha256']+'  '+x['path']+'\n' for x in files))
        write_json(d/'.transfer_identity.json',dict(transfer_id=tid,source_roots=manifest['source_roots'],manifest_sha256=sha(d/'transfer_manifest.json'),sums_sha256=sha(d/'SHA256SUMS')))
        (d/'rsync-exclude.txt').write_text('# Selection documentation. Explicit file whitelist overrides these patterns.\n# Local copying uses --files-from; remote upload copies the frozen selection.\n'+'\n'.join(sorted(EXCLUDE_NAMES)+[x+'/' for x in BUILD_DIR_PATTERNS]+PATTERNS)+
            '\n/mirheo_starter/runs/\n/mirheo_starter/results/\n/mirheo_starter/outputs/\n/mirheo_starter/logs/\n/mirheo_starter/test_code/outputs/\n/hemocell_starter/runs/\n/hemocell_starter/results/\n/hemocell_starter/outputs/\n/hemocell_starter/logs/\n')
        shutil.copytree(BASE/'tools',d/'metadata/transfer_tools',ignore=shutil.ignore_patterns('__pycache__'))
        write_json(d/'metadata/remote_verification.json',dict(status='PENDING'))
        (d/'metadata/transfer_log.txt').write_text(now()+' Local source snapshot prepared; source-to-copy hashes and execute bits matched.\n')
        write_json(d/'metadata/size_before_upload.json',packet_stats(d))
        state(d,'PREPARED',phase='FROZEN',source_regular_files=len(files),source_regular_bytes=manifest['regular_bytes'],symlinks=len(links),runtime_inputs_partial=True)
        say(json.dumps(dict(status='PREPARED',transfer_id=tid,upload_dir=str(d),statistics=packet_stats(d)),ensure_ascii=False))
    except Exception as e:
        state(d,'TRANSFER_INCOMPLETE',phase='LOCAL_PREPARATION_FAILED',error=str(e));raise

def validate_frozen(d):
    ident=read(d/'.transfer_identity.json')
    assert sha(d/'transfer_manifest.json')==ident['manifest_sha256'] and sha(d/'SHA256SUMS')==ident['sums_sha256'],'Frozen manifest changed'
    manifest=read(d/'transfer_manifest.json');expected={x['path']:x for x in manifest['files']}
    expected_links={x['path']:x for x in manifest['symlinks']};actual={};links={}
    for root,dirs,files in os.walk(d/'source',followlinks=False):
        for name in dirs+files:
            p=Path(root)/name;relative=p.relative_to(d/'source').as_posix()
            if p.is_symlink():links[relative]=p
            elif p.is_file():actual[relative]=p
            elif not p.is_dir():raise RuntimeError('Special file in frozen snapshot')
    assert set(actual)==set(expected) and set(links)==set(expected_links),'Frozen source file set changed'
    for key,p in actual.items():
        x=expected[key];s=p.stat()
        assert s.st_size==x['size'] and s.st_mode&0o111==x['execute_bits'] and sha(p)==x['sha256'],'Frozen source file mismatch: '+key
    for key,p in links.items():
        assert os.readlink(p)==expected_links[key]['target'] and not os.path.isabs(os.readlink(p)) and p.resolve(strict=True).is_relative_to(d/'source'),'Frozen symlink mismatch: '+key
    return manifest

def verify_remote(d):
    manifest=validate_frozen(d);paths=packet(d)
    before=read(d/'metadata/remote_before.json')['venv']
    return remote('verify',d,manifest_sha256=sha(d/'transfer_manifest.json'),sums_sha256=sha(d/'SHA256SUMS'),
        venv_before=before,critical_files=manifest['critical_files'],critical_directories=manifest['critical_directories'],
        packet_paths=paths,control_hashes={p:sha(d/p) for p in paths if not p.startswith('source/') and not (d/p).is_symlink()})

def sync(d,dry=False,control_only=False):
    argv=RSYNC+['-z','--partial','--stats','--info=progress2','--human-readable','-e',shlex.join(SSH)]
    if dry:argv+=['--dry-run']
    if control_only:
        listing=BASE/'runtime_logs'/(d.name+'-control-files.nul')
        listing.write_bytes(b''.join(p.encode()+b'\0' for p in packet(d) if not p.startswith('source/')))
        argv+=['--from0','--files-from='+str(listing)]
    argv += [str(d)+'/','vast-mirheo:/workspace/bloodflow/uploads/'+d.name+'/']
    label='dry-run' if dry else 'control-sync' if control_only else 'upload'
    logfile=BASE/'runtime_logs'/(d.name+'-'+label+'-'+datetime.now(timezone.utc).strftime('%H%M%S%f')+'.log')
    record=command_log(argv,logfile)
    say(label+' rsync exit_code='+str(record['exit_code'])+' log='+str(logfile))
    return record

def upload(d):
    local_check();current=read(d/'transfer_status.json')
    if current['status'] in VERIFIED:
        say('Completed snapshot: verifying without rewriting or regenerating it.');return verify(d)
    validate_frozen(d)
    try:inspection=remote('inspect')
    except Exception as e:state(d,'BLOCKED_BY_AUTHENTICATION',error=str(e));raise
    baseline=d/'metadata/remote_before.json'
    if baseline.exists():
        assert read(baseline)['venv']==inspection['venv'],'Existing remote .venv differs from initial transfer baseline; no overwrite performed'
    else:write_json(baseline,inspection)
    statistics=packet_stats(d)
    guard=dict(recorded_at=now(),**statistics,workspace_free_bytes=inspection['disk']['free'],
        projected_free_bytes=inspection['disk']['free']-statistics['regular_bytes'],staging_limit_bytes=3*1024**3,min_projected_free_bytes=5*1024**3)
    write_json(d/'metadata/size_before_upload.json',guard)
    say('UPLOAD_SIZE_GUARD '+json.dumps(guard,ensure_ascii=False))
    if guard['regular_bytes']>3*1024**3 or guard['projected_free_bytes']<5*1024**3:
        state(d,'BLOCKED_BY_SIZE_OR_DISK',phase='BEFORE_UPLOAD');return
    # Scan transfer-generated control records as well; source was scanned and then frozen.
    audit=scan_sensitive({p:d/p for p in packet(d) if not p.startswith('source/')})
    write_json(BASE/'runtime_logs'/(d.name+'-metadata-sensitive.json'),audit)
    if audit['status']!='PASS':raise RuntimeError('Control metadata sensitive-information audit requires review')
    remote('reserve',d)
    state(d,'TRANSFER_INCOMPLETE',phase='UPLOADING',last_error=None)
    dry=sync(d,dry=True)
    with (d/'metadata/transfer_log.txt').open('a') as f:f.write(json.dumps(dry,ensure_ascii=False)+'\n')
    if dry['exit_code']:state(d,'TRANSFER_INCOMPLETE',phase='DRY_RUN_FAILED');return
    say('Dry-run statistics:\n'+Path(dry['log_path']).read_text()[-2600:])
    actual=sync(d)
    with (d/'metadata/transfer_log.txt').open('a') as f:
        f.write(json.dumps(actual,ensure_ascii=False)+'\n'+Path(actual['log_path']).read_text()[-3000:]+'\n')
    if actual['exit_code']:state(d,'TRANSFER_INCOMPLETE',phase='RSYNC_FAILED');return
    # Deliver the completed rsync log before verifying the packet and its controls.
    controls=sync(d,control_only=True)
    if controls['exit_code']:state(d,'TRANSFER_INCOMPLETE',phase='LOG_SYNC_FAILED');return
    result=verify_remote(d)
    result.update(verified_at=now(),verification_scope='Full source SHA-256, counts, sizes, executable bits, links, exact packet file set, control hashes, existing venv fingerprint')
    write_json(d/'metadata/remote_verification.json',result)
    status='TRANSFER_VERIFIED_INPUTS_PARTIAL' if current.get('runtime_inputs_partial') else 'TRANSFER_VERIFIED'
    state(d,status,phase='REMOTE_VERIFIED',verified_at=result['verified_at'],sha256sum_exit_code=0,venv_unchanged=True,remote_free_bytes=result['disk']['free'])
    with (d/'metadata/transfer_log.txt').open('a') as f:f.write(now()+' UPLOAD_VERIFY_OK; real sha256sum exit=0; full source and packet checks passed.\n')
    finalsync=sync(d,control_only=True)
    if finalsync['exit_code']:state(d,'TRANSFER_INCOMPLETE',phase='FINAL_RECEIPT_SYNC_FAILED');return
    # Read-only final check also verifies delivered receipts; no recursive receipt rewriting.
    final=verify_remote(d)
    final.update(verified_at=now(),packet_statistics=packet_stats(d),final_receipts_synced=True)
    write_json(BASE/'runtime_logs'/(d.name+'-final-delivery.json'),final)
    say('UPLOAD_VERIFY_OK\n'+json.dumps(dict(transfer_id=d.name,status=status,regular_files=final['regular_files'],regular_bytes=final['regular_bytes'],symlinks=final['symlinks'],packet_statistics=final['packet_statistics'],disk=final['disk']),ensure_ascii=False))

def verify(d):
    local_check();result=verify_remote(d);result['verified_at']=now()
    logfile=BASE/'runtime_logs'/(d.name+'-reverify-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'.json')
    write_json(logfile,result);say('UPLOAD_VERIFY_OK\n'+json.dumps(dict(result=result,receipt=str(logfile)),ensure_ascii=False))

def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    sub.add_parser('prepare',help='Create a new frozen local snapshot; refuses if an unfinished transfer exists')
    for action in ['upload','verify','status']:
        p=sub.add_parser(action);p.add_argument('transfer_id',type=checked_id)
    args=parser.parse_args()
    if args.action=='prepare':prepare()
    else:
        d=directory(args.transfer_id)
        if args.action=='upload':upload(d)
        elif args.action=='verify':verify(d)
        else:say(json.dumps(read(d/'transfer_status.json'),ensure_ascii=False,indent=2))

if __name__=='__main__':
    try:main()
    except Exception as e:
        say('TRANSFER_STOPPED: '+str(e));sys.exit(2)

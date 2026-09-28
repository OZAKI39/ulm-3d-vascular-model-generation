"""Read-only source, result-state, remote mapping and publication checks."""
from pathlib import Path
import ast,csv,hashlib,json,os,re,subprocess,time
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260929';R=D/'brava_flow_roi_18mlmin'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    failures=[];changed_sources=[];trans={r['path']:r for r in read(M/'documentation_transforms.json')}
    rows=list(csv.DictReader((M/'source_files.csv').open()))
    for row in rows:
        source=Path(row['source']);target=D/row['destination'];h=row['sha256']
        if not source.exists() or sha(source)!=h:changed_sources.append(row['source'])
        if row['destination'] in trans:
            transform=trans[row['destination']]
            if sha(D/transform['original_copy'])!=h or sha(target)!=transform['snapshot_sha256']:failures.append('documentation transform: '+str(target))
        elif not target.is_file() or sha(target)!=h:failures.append('source copy: '+str(target))
    print('SOURCE_CHECK',len(rows),'changed sources',len(changed_sources),flush=True)
    server=read(M/'server_file_map.json')
    for row in server:
        if row['action'] not in {'MATCHED_LOCAL','COPIED_SERVER'}:continue
        if sha(D/row['destination'])!=row['sha256']:failures.append('server map: '+row['server_path'])
    active=read(R/'reports/ACTIVE_FLOW.json');assert active['status']=='PASS'
    case=R/active['case'];assert sha(case/'frozen_flow/steady_flow.vtu')==active['flow_sha256']
    for rel,h in read(case/'input_hashes.json').items():assert sha(case/rel)==h,rel
    for rel,h in active['postprocess_data_sha256'].items():assert sha(case/rel)==h,rel
    assert read(R/'gpu_solver_fix/backend_equivalence.json')['passed']
    assert read(R/'reports/LOCAL_FLOW_VERIFICATION.json')['PASS']
    videos=0
    for pose in ['0','15']:
        folder=R/'visualization'/('candidate_'+pose);media=read(folder/'MEDIA_VALIDATION.json')
        assert media['source_flow_sha256']==active['flow_sha256']
        for item in media['media']:
            assert item['frames']==432 and sha(folder/item['file'])==item['sha256'];videos+=1
    assert videos==10
    mb=R/'microbubble';stop=read(mb/'data/USER_STOP.json');assert stop['status']=='CANCELLED_BY_USER'
    assert not (mb/'data/TRAJECTORY_SEQUENCE_COMPLETE.json').exists()
    assert not (mb/'DELIVERY_COMPLETE.json').exists()
    assert not list((mb/'tracks').glob('mb_*/COMPLETE.json'))
    assert len(list((mb/'tracks').glob('mb_*')))==8
    assert not (mb/'data/query_cache_parity.json').exists()
    assert read(mb/'data/native_kernel_parity.json')['PASS']
    protected=read(mb/'data/preparation.json')['protected_sources']
    for rel,h in protected.items():assert sha(D/'ulm_particle_formal_p9a5'/rel)==h,rel
    pilot_files=0
    for kind in ['reference','native']:
        folder=mb/'pilot'/kind/'mb_000001'
        for rel,h in read(folder/'COMPLETE.json')['files'].items():assert sha(folder/rel)==h;pilot_files+=1
    # No protected scientific source is altered by this collection. Parse new
    # task orchestration only, without importing or starting its workflows.
    parsed=[]
    for folder in [R/'scripts',R/'microbubble/scripts',R/'visualization',M/'tools']:
        for p in folder.glob('*.py'):ast.parse(p.read_text());parsed.append(str(p.relative_to(D)))
    for row in read(M/'symlinks.json'):
        p=D/row['destination']
        if not p.is_symlink() or os.readlink(p)!=row['snapshot_target']:failures.append('link text: '+str(p))
        if not p.exists():failures.append('unresolved packaged link: '+str(p))
    patterns={
        'private_key':re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'),
        'github_token':re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,})'),
        'aws_access_key':re.compile(rb'AKIA[A-Z0-9]{16}'),
        'openai_key':re.compile(rb'sk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{40,}')}
    paths=set(subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=D).decode().split('\0'))
    paths.update(subprocess.check_output(['git','diff','--name-only','-z'],cwd=D).decode().split('\0'))
    findings=[];oversized=[]
    for rel in sorted(paths-{''}):
        p=D/rel
        if not p.is_file() or p.is_symlink():continue
        if p.stat().st_size>=100*1024**2:oversized.append(rel)
        if p.stat().st_size>8*1024**2 or p.suffix.lower() in {'.npz','.vtu','.vtp','.stl','.png','.pdf','.mp4','.gz','.zip'}:continue
        data=p.read_bytes()
        for name,pattern in patterns.items():
            if pattern.search(data):findings.append(dict(path=rel,pattern=name))
    result=dict(PASS=not failures and not changed_sources and not findings and not oversized,
        source_files_checked=len(rows),changed_sources=changed_sources,server_mappings_checked=sum(r['action'] in {'MATCHED_LOCAL','COPIED_SERVER'} for r in server),
        protected_scientific_sources_verified=len(protected),pilot_files_hash_verified=pilot_files,
        flow_sha256=active['flow_sha256'],flow_movies_hash_verified=videos,
        microbubble_status=stop['status'],production_completed_tracks=0,production_incomplete_tracks=8,
        query_cache_draft_validation='CANCELLED_NOT_PASSED',syntax_only_files=len(parsed),symlinks_checked=len(read(M/'symlinks.json')),
        credential_pattern_findings=findings,oversized_files=oversized,failures=failures,
        CFD_runs=0,new_trajectory_integrations=0,media_decode='Previously executed LOCAL_FLOW_VERIFICATION.json; this sync checks identical bytes, not new decoding',
        existing_mouse_validation='Retained from current_20260928; not rerun by this sync')
    (M/'validation/snapshot_validation.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,indent=2,ensure_ascii=False),flush=True)
    assert result['PASS']
if __name__=='__main__':main()

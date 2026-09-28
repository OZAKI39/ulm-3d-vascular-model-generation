"""Validate copied files, latest production identities, links and publication payload."""
from pathlib import Path
import collections,csv,hashlib,json,os,re,stat
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260928'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def save(name,value):(M/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def main():
    rows=list(csv.DictReader((M/'local_file_inventory.csv').open()))
    transformations={r['destination']:r for r in read(M/'documentation_transforms.json')}
    verified={};fail=[]
    for i,r in enumerate(rows):
        src=Path(r['source']);dst=D/r['destination'];want=r['sha256']
        if sha(src)!=want:fail.append(dict(path=str(src),reason='source changed during collection'))
        h=sha(dst);verified[str(dst.relative_to(D))]=h
        if h!=transformations.get(r['destination'],{}).get('snapshot_sha256',want):fail.append(dict(path=str(dst),reason='copy bytes differ'))
        if (i+1)%10000==0:print('source_and_snapshot_checked',i+1,flush=True)
    maps=read(M/'server_file_map.json');pointer_cache={};server_count=0
    for r in maps:
        dest=r['snapshot']
        if r['disposition']=='upstream_lfs_pointer_metadata_only':
            if dest not in pointer_cache:pointer_cache[dest]={x['sha256'] for x in read(D/dest)}
            assert r['sha256'] in pointer_cache[dest];continue
        if dest not in verified:verified[dest]=sha(D/dest)
        assert verified[dest]==r['sha256'],r['source'];server_count+=1
    links=read(M/'symlink_map.json')
    for r in links:
        p=D/r['destination'];assert p.is_symlink() and os.readlink(p)==r['snapshot_target']
        if not p.exists():fail.append(dict(path=r['destination'],reason='broken snapshot link'))
    project=D/'ulm_particle_formal_p9a5';run=project/'particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500'
    prod=read(run/'data/final_summary.json');delivery=read(run/'DELIVERY_COMPLETE.json')
    for name,want in delivery['files'].items():assert sha(run/name)==want,name
    for name,want in prod['identity']['protected_sources'].items():assert sha(project/name)==want,name
    statuses=collections.Counter();outlets=collections.Counter();track_files=0
    for folder in sorted((run/'tracks').glob('mb_*')):
        marker=read(folder/'COMPLETE.json');assert marker['identity']==prod['identity']
        for name,want in marker['files'].items():assert sha(folder/name)==want;track_files+=1
        metrics=read(folder/'metrics.json');statuses[metrics['status']]+=1
        if metrics['outlet']:outlets[metrics['outlet']]+=1
    assert sum(statuses.values())==1500 and statuses==prod['statuses'] and outlets==prod['outlets']
    stop=read(run/'stopping_analysis/data/analysis_summary.json');assert stop['all_pass'] and stop['total']==1500
    assert stop['status_counts']==dict(statuses) and stop['flow_sha256']==prod['flow_sha256']
    flow=D/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_ROI_only_balanced_pressure.vtu'
    assert sha(flow)==prod['flow_sha256']
    # Publication payload inventory. Do not recursively follow compatibility links.
    manifest=[];secret_findings=[];counts=collections.Counter();total=0;maxsize=0
    patterns=[('private_key',re.compile(rb'-----BEGIN (?:OPENSSH |RSA |EC |DSA )?PRIVATE KEY-----')),
              ('github_token',re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{60,})\b')),
              ('aws_access_id',re.compile(rb'\bAKIA[A-Z0-9]{16}\b')),
              ('openai_key',re.compile(rb'\bsk-(?:proj-)?[A-Za-z0-9_-]{45,}\b'))]
    ignored={str(M/'snapshot_manifest.csv'),str(M/'validation/snapshot_validation.json')}
    for root,dirs,files in os.walk(D,followlinks=False):
        dirs[:]=[n for n in dirs if n not in {'.git','__pycache__','.pytest_cache','.cache','mplconfig'}]
        for n in list(dirs):
            if (Path(root)/n).is_symlink():files.append(n);dirs.remove(n)
        for n in sorted(files):
            p=Path(root)/n
            if n=='.git' or n.endswith(('.pyc','.pyo')) or str(p) in ignored:continue
            rel=str(p.relative_to(D));islink=p.is_symlink()
            data=os.readlink(p).encode() if islink else None
            size=len(data) if islink else p.stat().st_size
            assert size<100*1024**2,rel
            h=hashlib.sha256();g=hashlib.sha1(b'blob '+str(size).encode()+b'\0')
            if islink:h.update(data);g.update(data)
            else:
                with p.open('rb') as f:
                    head=f.read(8192);textfile=b'\0' not in head
                    f.seek(0);tail=b''
                    while block:=f.read(1024**2):
                        h.update(block);g.update(block)
                        if textfile:
                            for label,pat in patterns:
                                if pat.search(tail+block):secret_findings.append(dict(path=rel,pattern=label))
                            tail=block[-256:]
            mode='120000' if islink else '100755' if p.stat().st_mode&stat.S_IXUSR else '100644'
            manifest.append(dict(path=rel,kind='symlink' if islink else 'file',bytes=size,sha256=h.hexdigest(),git_blob_oid=g.hexdigest(),mode=mode))
            total+=size;maxsize=max(maxsize,size);counts[rel.split('/')[0]]+=1
    with (M/'snapshot_manifest.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(manifest[0]));w.writeheader();w.writerows(manifest)
    result=dict(PASS=not fail and not secret_findings,source_files_verified=len(rows),server_content_mappings_verified=server_count,
        source_files_unchanged=True if not fail else False,symlinks_checked=len(links),track_files_verified=track_files,
        delivery_files_verified=len(delivery['files']),protected_source_files_verified=len(prod['identity']['protected_sources']),
        statuses=dict(statuses),outlets=dict(outlets),flow_sha256=prod['flow_sha256'],
        payload_files_and_links=len(manifest),payload_bytes=total,maximum_blob_bytes=maxsize,root_counts=dict(counts),
        secret_pattern_findings=secret_findings,failures=fail,CFD_runs=0,new_trajectory_integrations=0)
    save('validation/snapshot_validation.json',result);print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    assert result['PASS'],'See path-only findings in validation JSON; no secret values are logged'
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Compare every pre-existing artifact, scientific input, and old FEM entry."""
import gzip,json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,compare_inventory,git_state,sha256,write_json,now
R=ROOT/'reports/sv1_3h'
history=json.loads(gzip.decompress((R/'history_baseline.json.gz').read_bytes()));changes=[]
for key,expected in history['files'].items():
    p=ROOT/key
    if not p.exists() and not p.is_symlink():changes.append({'path':key,'reason':'missing'});continue
    s=p.lstat()
    if p.is_symlink():actual={'symlink':os.readlink(p),'mode':s.st_mode}
    elif p.is_file():actual={'sha256':sha256(p),'size':s.st_size,'mtime_ns':s.st_mtime_ns,'mode':s.st_mode}
    else:actual={'directory':True,'mode':s.st_mode}
    if actual!=expected:changes.append({'path':key,'before':expected,'after':actual})
for scope in history['scopes']:
    for key in inventory(ROOT/scope)['files']:
        if scope+'/'+key not in history['files']:changes.append({'path':scope+'/'+key,'reason':'added historical artifact'})
old=json.loads(gzip.decompress((R/'old_fem_baseline.json.gz').read_bytes()));current=inventory(old['root'])
old_result=compare_inventory(old,current);old_result['git_unchanged']=old['git']==git_state(old['root'])
if not old_result['git_unchanged']:old_result['status']='FAIL'
ref=json.loads((R/'reference_manifest.json').read_text())
ref_changed=[f['path'] for f in ref['files'] if sha256(ROOT/f['path'])!=f['sha256']]
official=git_state(ROOT/'external/svMultiPhysics')
official_result={'status':'PASS' if official==ref['official_source_git'] else 'FAIL','current':official}
d={'timestamp':now(),'status':'PASS' if not changes and not ref_changed and old_result['status']==official_result['status']=='PASS' else 'FAIL','historical':{'status':'FAIL' if changes else 'PASS','entries':len(history['files']),'changes':changes},'old_FEM':dict(old_result,entries=len(current['files'])),'reference_changed':ref_changed,'official_source':official_result,'atime_excluded':True}
write_json(R/'preservation_audit.json',d);print(json.dumps(d,ensure_ascii=False));raise SystemExit(0 if d['status']=='PASS' else 1)

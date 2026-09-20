#!/usr/bin/env python3
import gzip,json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import REPORT
from sv_validation.provenance import inventory,compare_inventory,git_state,sha256,write_json

history=json.loads(gzip.decompress((REPORT/'history_baseline.json.gz').read_bytes()))
checked={};changes=[]
for key,value in history['files'].items():
    p=ROOT/key
    if not p.exists() and not p.is_symlink():changes.append({'path':key,'reason':'missing'});continue
    stat=p.lstat()
    if p.is_symlink():actual={'symlink':os.readlink(p),'mode':stat.st_mode}
    elif p.is_file():actual={'sha256':sha256(p),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'mode':stat.st_mode}
    else:actual={'directory':True,'mode':stat.st_mode}
    checked[key]=actual
    if actual!=value:changes.append({'path':key,'before':value,'after':actual})
for folder in ('reports/sv1','outputs/sv1','logs/sv1','inputs'):
    for key in inventory(ROOT/folder)['files']:
        full=folder+'/'+key
        if full not in history['files']:changes.append({'path':full,'reason':'added historical artifact'})
old=json.loads(gzip.decompress((REPORT/'old_fem_baseline.json.gz').read_bytes()))
current=inventory(old['root']);comparison=compare_inventory(old,current)
result={'history_status':'PASS' if not changes else 'FAIL','historical_entries_checked':len(checked),'history_changes':changes,
        'old_fem_status':comparison['status'],'old_fem_entries_checked':len(current['files']),'old_fem_comparison':comparison,
        'old_fem_git_after':git_state(old['root'])}
if 'git' in old:result['old_fem_git_unchanged']=old['git']==result['old_fem_git_after']
result['status']='PASS' if not changes and comparison['status']=='PASS' and result.get('old_fem_git_unchanged',True) else 'FAIL'
write_json(REPORT/'preservation_audit.json',result);print(json.dumps(result,indent=2))
raise SystemExit(0 if result['status']=='PASS' else 1)

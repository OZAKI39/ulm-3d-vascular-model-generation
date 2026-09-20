"""Read-only historical, production, source and old-FEM preservation check."""
import gzip,json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,compare_inventory,git_state,sha256,now
R=ROOT/'reports/sv1_3p'
baseline=json.loads(gzip.decompress((R/'history_baseline.json.gz').read_bytes()))
inherited=json.loads(gzip.decompress((ROOT/'reports/sv1_3o/history_baseline.json.gz').read_bytes()))
expected=dict(inherited['files']);expected.update(baseline['files']);changes=[]
for n,e in expected.items():
 p=ROOT/n
 if not p.exists() and not p.is_symlink():changes.append(dict(path=n,reason='missing'));continue
 s=p.lstat()
 if p.is_symlink():actual=dict(symlink=os.readlink(p),mode=s.st_mode)
 elif p.is_file():actual=dict(sha256=sha256(p),size=s.st_size,mtime_ns=s.st_mtime_ns,mode=s.st_mode)
 else:actual=dict(directory=True,mode=s.st_mode)
 if actual!=e:changes.append(dict(path=n,reason='changed',before=e,after=actual))
for scope in set(baseline['scopes']+inherited['scopes']):
 for key in inventory(ROOT/scope)['files']:
  if '__pycache__' not in key and scope+'/'+key not in expected:changes.append(dict(path=scope+'/'+key,reason='added historical entry'))
old=json.loads(gzip.decompress((R/'old_fem_baseline.json.gz').read_bytes()));current=inventory(old['root']);old_result=compare_inventory(old,current)
old_result['git_unchanged']=git_state(old['root'])==old['git']
ref=json.loads((R/'reference_freeze.json').read_text());official=git_state(ROOT/'external/svMultiPhysics')
source=json.loads((R/'source_patch.json').read_text());source_changes=[n for n,h in source['after'].items() if sha256(ROOT/source['source']/n)!=h]
d=dict(status='PASS' if not changes and old_result['status']=='PASS' and old_result['git_unchanged'] and official==ref['official_source_git'] and not source_changes else 'FAIL',timestamp=now(),historical_entries=len(expected),historical_changes=changes,old_FEM=old_result,official_source_unchanged=official==ref['official_source_git'],adopted_source_unexpected_changes=source_changes,production='CPU_EARLY_STOP_PRODUCTION',atime_excluded=True)
(R/'preservation_audit.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d));raise SystemExit(0 if d['status']=='PASS' else 1)

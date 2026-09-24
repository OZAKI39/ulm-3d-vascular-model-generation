from pathlib import Path
import hashlib,json,subprocess,sys,time
from concurrent.futures import ThreadPoolExecutor
root=Path(sys.argv[1]);start=time.time()
manifest=json.loads((root/'data/taylor_hood_remote_protection_before.json').read_text())
def check(item):
 p,before=item;p=Path(p)
 if not p.is_file():return dict(path=str(p),reason='MISSING')
 with p.open('rb') as f:after=hashlib.file_digest(f,'sha256').hexdigest()
 if after!=before['sha256'] or p.stat().st_size!=before['size']:return dict(path=str(p),before=before,after_sha256=after)
with ThreadPoolExecutor(max_workers=6) as pool:changed=[v for v in pool.map(check,manifest['files'].items()) if v]
git={}
for args in [('status','--short'),('diff','--stat'),('rev-parse','HEAD')]:
 p=subprocess.run(['git','--no-optional-locks','-C',manifest['source_path'],*args],capture_output=True,text=True)
 git[' '.join(args)]=dict(stdout=p.stdout,stderr=p.stderr,returncode=p.returncode)
result=dict(status='PASS' if not changed else 'CONCURRENT_WORKTREE_CHANGE',files_checked=len(manifest['files']),files_changed=len(changed),changes=changed,git_after=git,restoration_attempted=False,elapsed_seconds=time.time()-start)
(root/'logs/taylor_hood_remote_protection_final.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k not in ['git_after','changes']}))

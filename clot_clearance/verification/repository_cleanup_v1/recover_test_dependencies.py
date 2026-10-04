"""Restore accidentally removed indirect manifest dependencies by exact hash."""
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import shutil
from audit_inventory import AUDIT, PROJECT, ROOT, dump_csv, dump_json


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    missing=json.loads((AUDIT/'TEST_MANIFEST_DEPENDENCIES_MISSING.json').read_text())
    rows=json.loads((AUDIT/'CLASSIFIED_INVENTORY.json').read_text());by_path={r['path']:r for r in rows}
    available=defaultdict(list)
    for r in rows:
        if r['type']=='file' and (ROOT/r['path']).is_file():available[r['hash']].append(ROOT/r['path'])
    for path,r in missing.items():
        assert by_path[path]['hash']==r['hash']
    needed={r['hash'] for r in missing.values() if not available[r['hash']]}
    sizes={by_path[p]['size_bytes'] for p,r in missing.items() if r['hash'] in needed}
    # Read-only search of local sibling workspaces and task/pytest scratch.
    candidates=[];errors=[]
    for base in [Path('/home/lzy/projects'),Path('/tmp')]:
        for directory,dirs,files in os.walk(base,followlinks=False,onerror=lambda e:errors.append(str(e))):
            dirs[:]=[d for d in dirs if d not in ('.git','node_modules') and Path(directory)/d!=ROOT]
            for name in files:
                p=Path(directory)/name
                try:
                    if not p.is_symlink() and p.stat().st_size in sizes and p.suffix in ('.json','.npz','.gz','.nbc','.nbi'):
                        candidates.append(p)
                except OSError:pass
    print('Searching',len(candidates),'same-size recovery candidates for',len(needed),'hashes',flush=True)
    def check(p):
        try:return p,digest(p)
        except OSError:return p,''
    with ThreadPoolExecutor(max_workers=4) as pool:
        for p,h in pool.map(check,candidates):
            if h in needed:available[h].append(p)
    plans=[];unavailable=[]
    for name,r in missing.items():
        if os.path.lexists(ROOT/name):
            assert digest(ROOT/name)==r['hash'];continue
        sources=available.get(r['hash'],[])
        if not sources:
            unavailable.append(dict(path=name,hash=r['hash'],size_bytes=by_path[name]['size_bytes'],dependency=r['dependency']));continue
        source=sorted(sources,key=lambda p:len(str(p)))[0]
        assert digest(source)==r['hash']
        plans.append(dict(path=name,source=str(source),size_bytes=by_path[name]['size_bytes'],hash=r['hash'],reason='Indirect existing test dependency on historical manifest; restore byte-identical canonical copy.'))
    dump_csv('TEST_DEPENDENCY_RESTORATION_MANIFEST.csv',plans,['path','source','size_bytes','hash','reason'])
    log=[]
    for r in plans:
        target=ROOT/r['path'];target.parent.mkdir(parents=True,exist_ok=True)
        with Path(r['source']).open('rb') as src,target.open('xb') as dst:shutil.copyfileobj(src,dst)
        assert digest(target)==r['hash']
        original=by_path[r['path']];os.utime(target,ns=(original['mtime_ns'],original['mtime_ns']))
        log.append(dict(**r,status='RESTORED_EXACT'))
    dump_csv('TEST_DEPENDENCY_RESTORATION_LOG.csv',log,['path','source','size_bytes','hash','reason','status'])
    dump_json('RECOVERY_UNAVAILABLE.json',unavailable)
    dump_json('RECOVERY_SEARCH.json',dict(candidates=len(candidates),needed_hashes=len(needed),restored=len(log),remaining=len(unavailable),read_only_search_errors=errors))
    print('Restored',len(log),'remaining',len(unavailable),flush=True)


if __name__=='__main__':main()

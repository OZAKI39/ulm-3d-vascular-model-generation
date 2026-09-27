"""Content provenance without importing or writing into the old FEM project."""
import hashlib
import json
import os
import subprocess
from datetime import datetime,timezone
from pathlib import Path

def now():return datetime.now(timezone.utc).isoformat()
def sha256(path):
    d=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):d.update(block)
    return d.hexdigest()
def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def inventory(root):
    root=Path(root).resolve();files={}
    for base,dirs,names in os.walk(root,followlinks=False):
        for name in sorted(dirs+names):
            p=Path(base)/name;s=p.lstat();key=str(p.relative_to(root))
            if p.is_symlink():files[key]={'symlink':os.readlink(p),'mode':s.st_mode}
            elif p.is_file():files[key]={'sha256':sha256(p),'size':s.st_size,'mtime_ns':s.st_mtime_ns,'mode':s.st_mode}
            elif p.is_dir():files[key]={'directory':True,'mode':s.st_mode}
    return {'root':str(root),'timestamp':now(),'files':files,'atime_excluded':True,'symlinks_followed':False}
def git_state(root):
    env=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
    def read(*args):return subprocess.check_output(['git','-C',str(root),*args],env=env,text=True)
    return {'head':read('rev-parse','HEAD').strip(),'status':read('status','--porcelain=v1','--untracked-files=all'),
            'diff_sha256':hashlib.sha256(read('diff','--binary').encode()).hexdigest(),
            'cached_diff_sha256':hashlib.sha256(read('diff','--cached','--binary').encode()).hexdigest()}
def compare_inventory(before,after):
    a=before['files'];b=after['files']
    changes={'modified':sorted(k for k in a.keys()&b.keys() if a[k]!=b[k]),'added':sorted(b.keys()-a.keys()),'deleted':sorted(a.keys()-b.keys())}
    return {'status':'FAIL' if any(changes.values()) else 'PASS',**changes}

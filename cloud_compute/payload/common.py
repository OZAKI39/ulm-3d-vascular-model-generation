"""Standard-library primitives shared by local orchestration and remote jobs."""
import hashlib,json,os,stat,time
from pathlib import Path,PurePosixPath
from datetime import datetime,timezone

def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+'.tmp')
    q.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n');q.replace(p)
def object_hash(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def safe_relative(value):
    if not isinstance(value,str) or not value or '\\' in value or '\x00' in value:raise ValueError('UNSAFE_MANIFEST_PATH')
    p=PurePosixPath(value)
    if p.is_absolute() or any(x in ('','.','..') for x in value.split('/')):raise ValueError('UNSAFE_MANIFEST_PATH')
    return value
def files(root):
    rows=[];root=Path(root)
    for parent,dirs,names in os.walk(root,followlinks=False):
        for name in sorted(dirs+names):
            p=Path(parent)/name;s=p.lstat();relative=p.relative_to(root).as_posix();safe_relative(relative)
            if p.is_symlink():raise ValueError('SYMLINK_NOT_ALLOWED_IN_RESULTS:'+relative)
            if stat.S_ISREG(s.st_mode):rows.append(dict(path=relative,size=s.st_size,sha256=sha(p),execute_bits=s.st_mode&0o111))
            elif not stat.S_ISDIR(s.st_mode):raise ValueError('SPECIAL_FILE_NOT_ALLOWED:'+relative)
    return sorted(rows,key=lambda x:x['path'])
def validate_result(root):
    root=Path(root)
    if root.is_symlink():raise ValueError('UNSAFE_RESULT_ROOT')
    ready=root/'RESULTS_READY';manifest=root/'results_manifest.json'
    if not ready.is_file() or not manifest.is_file():raise ValueError('RESULTS_NOT_READY')
    mark=read(ready);m=read(manifest)
    if mark['manifest_sha256']!=sha(manifest):raise ValueError('MANIFEST_HASH_MISMATCH')
    if not m.get('files_closed') or m.get('state') not in ('COMPLETED','FAILED','BLOCKED','TIMEOUT','DISK_LIMIT'):raise ValueError('RESULTS_NOT_CLOSED')
    expected={}
    for row in m['files']:
        p=safe_relative(row['path'])
        if p in expected or row.get('type','file')!='file':raise ValueError('UNSAFE_MANIFEST_ENTRY')
        expected[p]=row
    actual={x['path']:x for x in files(root) if x['path'] not in ('RESULTS_READY','results_manifest.json','LOCAL_ARCHIVE_VERIFIED.json')}
    if set(actual)!=set(expected):raise ValueError('RESULT_FILE_SET_MISMATCH')
    for p,row in actual.items():
        if any(row[k]!=expected[p][k] for k in ('size','sha256','execute_bits')):raise ValueError('RESULT_CONTENT_MISMATCH:'+p)
    return dict(status='RESULTS_RETURN_VERIFIED',job_id=m['job_id'],run_state=m['state'],files=len(actual),bytes=sum(x['size'] for x in actual.values()),manifest_sha256=sha(manifest))
def seal(root,job_id,state):
    root=Path(root)
    assert not (root/'RESULTS_READY').exists(),'Never reseal completed results'
    write(root/'results_manifest.json',dict(job_id=job_id,state=state,files_closed=True,created_at=now(),files=files(root)))
    write(root/'RESULTS_READY',dict(job_id=job_id,manifest_sha256=sha(root/'results_manifest.json'),ready_at=now()))
def disk_allowed(free,kind,growth=0):
    floor=5*1024**3 if kind=='build' else 3*1024**3
    return free-growth>=floor

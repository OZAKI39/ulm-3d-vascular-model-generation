"""Standard-library primitives shared by local orchestration and remote jobs."""
import hashlib,json,os,stat,time,signal,subprocess
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
def tree_identity(root):
    """Read-only content/mode/mtime identity; preserve and never traverse symlinks."""
    root=Path(root);rows=[]
    for parent,dirs,names in os.walk(root,followlinks=False):
        for name in sorted(dirs+names):
            p=Path(parent)/name;s=p.lstat();row=dict(path=p.relative_to(root).as_posix(),mode=s.st_mode,mtime_ns=s.st_mtime_ns)
            if p.is_symlink():row['link']=os.readlink(p)
            elif p.is_file():row.update(size=s.st_size,sha256=sha(p))
            rows.append(row)
    return dict(files=len(rows),fingerprint=object_hash(sorted(rows,key=lambda x:x['path'])))
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
    if not m.get('files_closed') or m.get('state') not in ('COMPLETED','FAILED','BLOCKED','TIMEOUT','DISK_LIMIT','STOPPED'):raise ValueError('RESULTS_NOT_CLOSED')
    expected={}
    for row in m['files']:
        p=safe_relative(row['path'])
        if p in expected or row.get('type','file')!='file':raise ValueError('UNSAFE_MANIFEST_ENTRY')
        expected[p]=row
    actual={x['path']:x for x in files(root) if x['path'] not in ('RESULTS_READY','results_manifest.json','LOCAL_ARCHIVE_VERIFIED.json')}
    if set(actual)!=set(expected):raise ValueError('RESULT_FILE_SET_MISMATCH')
    for p,row in actual.items():
        if any(row[k]!=expected[p][k] for k in ('size','sha256','execute_bits')):raise ValueError('RESULT_CONTENT_MISMATCH:'+p)
    if (root/'run_status.json').is_file() and read(root/'run_status.json').get('protocol')=='CLOUD_FULL_RBC_V1':
        state=read(root/'run_status.json')
        for name in ['actual_parameters.json','provenance.json','execution.json','full_spec.json']:
            if not (root/name).is_file():raise ValueError('MISSING_FULL_RBC_RECORD:'+name)
            read(root/name)
        if state['status']=='COMPLETED':
            c=read(root/'simulation/completion.json');s=read(root/'full_spec.json')
            if not c.get('completed') or c.get('shear_steps')!=400000 or c.get('prep_steps')!=60000 or state.get('exit_code')!=0 or abs(c.get('strain_end',-1)-4)>1e-5:
                raise ValueError('FALSE_FULL_RBC_COMPLETION')
            if s.get('cloud_full_protocol')!='CLOUD_FULL_RBC_V1' or s.get('dt')!=.0005:raise ValueError('WRONG_FULL_RBC_PROTOCOL')
            import csv,math
            with (root/'simulation/vertices.csv').open() as stream:rows=[r for r in csv.DictReader(stream) if r['phase']=='shear' and int(r['step'])==460000]
            if len(rows)!=642 or not all(math.isfinite(float(r[k])) for r in rows for k in ['x','y','z']):raise ValueError('FULL_RBC_ENDPOINT_DATA_MISSING')
    if (root/'run_status.json').is_file() and read(root/'run_status.json').get('status')=='CLOUD_SMOKE_PASS':
        import csv,math
        c=read(root/'B_liquid/completion.json')
        if c.get('actual_steps',0)<=0 or not c.get('finite'):raise ValueError('EMPTY_OR_INVALID_LIQUID_RESULT')
        with (root/'B_liquid/particles_final.csv').open() as f:rows=list(csv.DictReader(f))
        if len(rows)!=c['particles'] or not rows:raise ValueError('EMPTY_PARTICLE_CSV')
        if not all(math.isfinite(float(row[k])) for row in rows for k in ['x','y','z','vx','vy','vz']):raise ValueError('NONFINITE_PARTICLE_CSV')
        if (root/'C_repair/completion.json').is_file():
            c=read(root/'C_repair/completion.json')
            if not c['completed'] or c['actual_steps']<=0 or c['membrane_updates']<=0:raise ValueError('INVALID_REPAIR_COMPLETION')
            with (root/'C_repair/vertices.csv').open() as f:rows=list(csv.DictReader(f))
            if len(rows)<c['vertices']:raise ValueError('EMPTY_MEMBRANE_CSV')
            if not all(math.isfinite(float(row[k])) for row in rows for k in ['x','y','z','fx','fy','fz']):raise ValueError('NONFINITE_MEMBRANE_CSV')
    return dict(status='RESULTS_RETURN_VERIFIED',job_id=m['job_id'],run_state=m['state'],files=len(actual),bytes=sum(x['size'] for x in actual.values()),manifest_sha256=sha(manifest))
def seal(root,job_id,state):
    root=Path(root)
    assert not (root/'RESULTS_READY').exists(),'Never reseal completed results'
    if state not in ('COMPLETED','FAILED','BLOCKED','TIMEOUT','DISK_LIMIT','STOPPED'):raise ValueError('CANNOT_SEAL_RUNNING_RESULTS')
    write(root/'results_manifest.json',dict(job_id=job_id,state=state,files_closed=True,created_at=now(),files=files(root)))
    write(root/'RESULTS_READY',dict(job_id=job_id,manifest_sha256=sha(root/'results_manifest.json'),ready_at=now()))
def disk_allowed(free,kind,growth=0):
    floor=5*1024**3 if kind=='build' else 3*1024**3
    return free-growth>=floor

def wait_guarded(process,seconds,free_bytes,grace=2,stop_reason=None):
    """Stop only the caller-owned process group; reserve termination grace within the limit."""
    start=time.monotonic();grace=min(grace,max(.01,seconds/4));reason=None
    while process.poll() is None:
        elapsed=time.monotonic()-start
        if elapsed>=seconds-grace:reason='TIMEOUT'
        elif not disk_allowed(free_bytes(),'run'):reason='DISK_LIMIT'
        elif stop_reason is not None:reason=stop_reason()
        if reason:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=max(.01,seconds-(time.monotonic()-start)))
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
            break
        time.sleep(min(.1,max(.005,seconds/20)))
    return dict(exit_code=process.returncode,elapsed_s=time.monotonic()-start,stop_reason=reason)

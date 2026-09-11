"""Run over SSH stdin using only system Python's standard library."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()


def tree(p):
    rows=[]
    for base,dirs,files in os.walk(p,followlinks=False):
        for name in sorted(dirs+files):
            q=Path(base)/name;s=q.lstat();relative=q.relative_to(p).as_posix()
            if stat.S_ISLNK(s.st_mode):rows.append([relative,'symlink',os.readlink(q)])
            elif stat.S_ISREG(s.st_mode):rows.append([relative,'file',s.st_size,stat.S_IMODE(s.st_mode),s.st_mtime_ns,digest(q)])
            elif stat.S_ISDIR(s.st_mode):rows.append([relative,'directory',stat.S_IMODE(s.st_mode)])
            else:rows.append([relative,'SPECIAL'])
    encoded=json.dumps(sorted(rows),separators=(',',':')).encode()
    return dict(path=str(p),entries=len(rows),regular_files=sum(x[1]=='file' for x in rows),
                content_and_mode_fingerprint=hashlib.sha256(encoded).hexdigest())


def main():
    request=json.loads(sys.argv[1]);mode=request['mode']
    workspace=Path('/workspace/bloodflow');venv=workspace/'.venv'
    assert workspace.is_dir() and venv.is_dir(),'EXPECTED_WORKSPACE_OR_VENV_MISSING'
    assert workspace.resolve()==workspace,'WORKSPACE_SYMLINK_NOT_ALLOWED'
    disk=shutil.disk_usage('/workspace')
    result=dict(mode=mode,user=__import__('pwd').getpwuid(os.getuid()).pw_name,uid=os.getuid(),hostname=__import__('socket').gethostname(),
                disk=dict(total=disk.total,used=disk.used,free=disk.free),rsync=shutil.which('rsync'),venv=tree(venv))
    if mode=='inspect':
        guides=[]
        for name in ['/AGENTS.md','/workspace/AGENTS.md','/workspace/bloodflow/AGENTS.md','/workspace/bloodflow/uploads/AGENTS.md']:
            p=Path(name)
            if p.is_file():guides.append(dict(path=str(p),sha256=digest(p)))
        result.update(status='PASS',instruction_files=guides);return result
    target=Path(request['remote_dir']);base=workspace/'uploads'
    assert target.parent==base and target.name==request['transfer_id'],'UNSAFE_TARGET'
    assert not base.is_symlink() and target.resolve()==target,'TARGET_SYMLINK_NOT_ALLOWED'
    identity=request['identity']
    if mode=='reserve':
        if target.exists():
            p=target/'.transfer_identity.json'
            assert p.is_file() and json.loads(p.read_text())==identity,'EXISTING_TARGET_IDENTITY_MISMATCH'
        else:
            target.mkdir(parents=True,exist_ok=False)
            (target/'.transfer_identity.json').write_text(json.dumps(identity,ensure_ascii=False,indent=2)+'\n')
        result.update(status='PASS',remote_dir=str(target));return result
    assert mode=='verify' and target.is_dir()
    assert json.loads((target/'.transfer_identity.json').read_text())==identity
    assert digest(target/'transfer_manifest.json')==request['manifest_sha256'],'REMOTE_MANIFEST_CHANGED'
    assert digest(target/'SHA256SUMS')==request['sums_sha256'],'REMOTE_CHECKSUM_LIST_CHANGED'
    manifest=json.loads((target/'transfer_manifest.json').read_text());source=target/'source'
    expected={x['path']:x for x in manifest['files']};expected_links={x['path']:x for x in manifest['symlinks']}
    actual={};links={};errors=[]
    for root,dirs,files in os.walk(source,followlinks=False):
        for name in dirs+files:
            p=Path(root)/name;s=p.lstat();relative=p.relative_to(source).as_posix()
            if stat.S_ISLNK(s.st_mode):links[relative]=p
            elif stat.S_ISREG(s.st_mode):actual[relative]=p
            elif not stat.S_ISDIR(s.st_mode):errors.append('SPECIAL_FILE:'+relative)
    if set(actual)!=set(expected):errors.append('REGULAR_FILE_SET_MISMATCH')
    if set(links)!=set(expected_links):errors.append('SYMLINK_SET_MISMATCH')
    for relative,p in actual.items():
        if relative not in expected:continue
        entry=expected[relative];s=p.stat()
        if s.st_size!=entry['size']:errors.append('SIZE:'+relative)
        if s.st_mode&0o111!=entry['execute_bits']:errors.append('EXECUTE_BITS:'+relative)
    for relative,p in links.items():
        if relative not in expected_links:continue
        value=os.readlink(p)
        if value!=expected_links[relative]['target'] or os.path.isabs(value):errors.append('LINK_TARGET:'+relative)
        try:
            dest=p.resolve(strict=True)
            if not dest.is_relative_to(source):errors.append('EXTERNAL_LINK:'+relative)
        except (OSError,RuntimeError):errors.append('BROKEN_OR_LOOPING_LINK:'+relative)
    for relative in request['critical_files']:
        p=source/relative
        if not p.is_file() or p.stat().st_size==0:errors.append('MISSING_REQUIRED_FILE:'+relative)
    for relative in request['critical_directories']:
        if not (source/relative).is_dir():errors.append('MISSING_REQUIRED_DIRECTORY:'+relative)
    packet=set()
    for root,dirs,files in os.walk(target,followlinks=False):
        for name in dirs+files:
            p=Path(root)/name
            if p.is_file() or p.is_symlink():packet.add(p.relative_to(target).as_posix())
    extra=sorted(packet-set(request['packet_paths']))
    missing=sorted(set(request['packet_paths'])-packet)
    if extra:errors.append('UNEXPECTED_PACKET_FILES')
    if missing:errors.append('MISSING_PACKET_FILES')
    for relative,expected_hash in request.get('control_hashes',{}).items():
        p=target/relative
        if not p.is_file() or p.is_symlink() or digest(p)!=expected_hash:errors.append('CONTROL_FILE_HASH:'+relative)
    if result['venv']!=request['venv_before']:errors.append('EXISTING_VENV_CHANGED')
    check=subprocess.run(['sha256sum','--quiet','-c','../SHA256SUMS'],cwd=source,capture_output=True,text=True)
    if check.returncode:errors.append('SHA256SUMS_FAILED')
    result.update(status='PASS' if not errors else 'FAIL',errors=errors,extra_packet_files=extra,missing_packet_files=missing,
        sha256sum_exit_code=check.returncode,sha256sum_diagnostics=(check.stdout+check.stderr)[:2000],
        regular_files=len(actual),regular_bytes=sum(p.stat().st_size for p in actual.values()),symlinks=len(links),
        execute_bits_verified=not any(x.startswith('EXECUTE_BITS:') for x in errors),
        venv_unchanged=result['venv']==request['venv_before'],remote_dir=str(target))
    if not errors and check.returncode==0:result['marker']='UPLOAD_VERIFY_OK'
    return result


if __name__=='__main__':
    try:
        result=main();print(json.dumps(result,ensure_ascii=False));sys.exit(0 if result['status']=='PASS' else 2)
    except Exception as e:
        print(json.dumps(dict(status='FAIL',error=type(e).__name__+': '+str(e)),ensure_ascii=False));sys.exit(2)

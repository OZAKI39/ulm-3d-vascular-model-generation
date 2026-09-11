"""Read-only selection and audit helpers for the WSL source upload."""
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tarfile
import zipfile

BASE=Path('/home/lzy/projects/cloud_upload')
PROJECTS=Path('/home/lzy/projects')
ROOTS={name:PROJECTS/name for name in ('mirheo_starter','hemocell_starter')}
CONFIG=json.loads((BASE/'tools/selection.json').read_text())
CAMPAIGN=CONFIG['repair_campaign']
B=ROOTS['mirheo_starter']/'data/single_rbc_repair'/CAMPAIGN
RUN=ROOTS['mirheo_starter']/'runs/single_rbc_repair'/CAMPAIGN/'gpu'/CONFIG['latest_run']
EXCLUDE_NAMES={'.git','.venv','venv','env','conda','miniconda3','anaconda3','.conda',
    '__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','.cache','analysis_cache',
    'CMakeFiles','downloads','download_cache','install_cache','.ssh','.codex','node_modules',
    '.tox','.nox','.ccache','.nv','_skbuild','dist'}
BUILD_DIR_PATTERNS=['build','build-*','build_*','cmake-build-*']
PATTERNS=['*.pyc','*.pyo','*.o','*.obj','*.a',
    '*.so','*.so.*','*.dll','*.exe','*.dylib','*.whl','*.deb','*.rpm','.env','.env.*',
    'id_rsa*','id_ed25519*','*.key','*.p12','*.pfx','.netrc','.npmrc','.pypirc',
    'credentials','credentials.json','*.swp','*.swo']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def write_json(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.tmp')
    tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n');tmp.replace(p)

def identity(p):
    s=p.lstat()
    return dict(device=s.st_dev,inode=s.st_ino,size=s.st_size,mtime_ns=s.st_mtime_ns,
                ctime_ns=s.st_ctime_ns,mode=stat.S_IMODE(s.st_mode))

def git_info(p):
    env=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
    def git(*args):
        r=subprocess.run(['git','-C',str(p),*args],capture_output=True,text=True,env=env)
        return r.stdout.strip() if r.returncode==0 else None
    top=git('rev-parse','--show-toplevel')
    exact=top is not None and Path(top).resolve()==p.resolve()
    return dict(path=str(p),is_exact_git_root=exact,detected_top_level=top,
        head=git('rev-parse','HEAD') if exact else None,
        status_porcelain=git('status','--porcelain=v1','--untracked-files=all') if exact else None,
        submodule_status=git('submodule','status','--recursive') if exact else None)

def whitelist():
    paths={p:'Build/source evidence exception' for p in CONFIG['supplementary_build_evidence']+CONFIG['upstream_source_under_build_whitelist']}
    old=ROOTS['mirheo_starter']/'runs/single_rbc_benchmark'/CONFIG['legacy_campaign']/'solver/main_hemocell_1'
    for name in CONFIG['case_whitelist']:
        p=old/name
        if p.is_file():paths[p.relative_to(PROJECTS).as_posix()]='Historical HemoCell input example and provenance; not a new cloud parameter freeze'
    for run in sorted(RUN.parent.iterdir()):
        if not run.is_dir():continue
        for name in ['spec.json','execution.json','completion.json','case_creation.json','worker_source.py','continuous_protocol_source.py','handoff.json']:
            p=run/name
            if p.is_file():paths[p.relative_to(PROJECTS).as_posix()]='Explicit A0-A6 source/library and case provenance; excludes raw result tree'
    for p in [RUN.parent/'budget_ledger.json',RUN.parent.parent/'build/budget_ledger.json']:
        if p.is_file():paths[p.relative_to(PROJECTS).as_posix()]='Historical local budget evidence only; not cloud execution authorization'
    return paths

def reason_for(project,relative):
    for item in CONFIG['optional_exclusions']:
        if project==item['project'] and (relative==item['path'] or relative.startswith(item['path']+'/')):
            return item['reason']
    parts=Path(relative).parts
    if parts[0] in ('runs','results','outputs','logs'):return 'Default root history exclusion; explicit file whitelist takes precedence'
    if relative=='test_code/outputs' or relative.startswith('test_code/outputs/'):return 'Generated test outputs'
    for index,name in enumerate(parts):
        if name in EXCLUDE_NAMES or any(fnmatch.fnmatchcase(name,p) for p in PATTERNS):return 'Default excluded environment/cache/build/binary/credential path: '+name
        if any(fnmatch.fnmatchcase(name,p) for p in BUILD_DIR_PATTERNS) and (index<len(parts)-1 or (ROOTS[project]/relative).is_dir()):return 'Default build directory exclusion: '+name
    return None

def inventory():
    included={};excluded=[];special=[];supplement=whitelist()
    def omit(p,key,reason,detailed=False):
        s=p.lstat();row=dict(path=key,source_path=str(p),reason=reason,type='directory' if p.is_dir() and not p.is_symlink() else 'symlink' if p.is_symlink() else 'file',size=s.st_size)
        if p.is_symlink():row.update(target=os.readlink(p),link_text_sha256=hashlib.sha256(os.readlink(p).encode()).hexdigest())
        elif detailed and p.is_file():row['sha256']=sha(p)
        excluded.append(row)
    for project,root in ROOTS.items():
        for parent,dirs,files in os.walk(root,followlinks=False):
            for name in sorted(dirs+files):
                p=Path(parent)/name;relative=p.relative_to(root).as_posix();key=project+'/'+relative
                why=reason_for(project,relative)
                if why:
                    detailed='analysis_cache' in Path(relative).parts or any(relative==i['path'] for i in CONFIG['optional_exclusions'] if i['project']==project)
                    omit(p,key,why,detailed)
                    if p.is_dir() and not p.is_symlink():
                        if detailed:
                            for dp,ds,fs in os.walk(p,followlinks=False):
                                for n in ds+fs:
                                    q=Path(dp)/n
                                    if not q.is_dir() or q.is_symlink():omit(q,q.relative_to(PROJECTS).as_posix(),why,True)
                        dirs.remove(name)
                    continue
                s=p.lstat()
                if stat.S_ISDIR(s.st_mode):continue
                if not (stat.S_ISREG(s.st_mode) or stat.S_ISLNK(s.st_mode)):
                    special.append(dict(path=key,type='SPECIAL'));continue
                if stat.S_ISREG(s.st_mode):
                    with p.open('rb') as f:magic=f.read(8)
                    if magic.startswith((b'\x7fELF',b'MZ',b'!<arch>\n',b'\xcf\xfa\xed\xfe',b'\xfe\xed\xfa\xcf')):
                        omit(p,key,'Compiled executable/library detected by magic, irrespective of execute bit',True);continue
                included[key]=dict(source_path=str(p),explicit_supplemental=False)
    for key,reason in supplement.items():
        p=PROJECTS/key
        if not p.is_file() or p.is_symlink():raise RuntimeError('Missing or unsafe explicit supplement: '+key)
        included[key]=dict(source_path=str(p),explicit_supplemental=True,supplement_reason=reason)
    return included,excluded,special

SECRET_RULES={
    'private_key_block':re.compile(rb'-----BEGIN (?:[A-Z0-9]+ )?PRIVATE KEY-----'),
    'provider_token':re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|sk-(?:proj-)?[A-Za-z0-9_-]{32,}|AKIA[A-Z0-9]{16}|xox[baprs]-[A-Za-z0-9-]{20,})\b'),
    'credential_url':re.compile(rb'https?://[^\s/<>"\x27:@]{1,80}:[^\s/<>"\x27@]{1,120}@'),
    'authorization_bearer':re.compile(rb'(?i)authorization["\x27\s:=]+bearer\s+[A-Za-z0-9._-]{16,}'),
    'credential_assignment':re.compile(rb'(?i)["\x27]?(?:api[_-]?key|access[_-]?token|auth[_-]?token|password|passwd|secret[_-]?key)["\x27]?\s*[:=]\s*["\x27]([^"\x27\r\n]{8,})["\x27]')}

def scan_sensitive(paths):
    findings=[];count=0;archive_members=0;archive_errors=[];file_hashes={}
    def scan(data,label):
        nonlocal count
        count+=1
        lower=data.lower()
        probes={'private_key_block':[b'private key'], 'provider_token':[b'ghp_',b'gho_',b'ghu_',b'ghs_',b'ghr_',b'github_pat_',b'sk-',b'akia',b'xox'],
            'credential_url':[b'http://',b'https://'], 'authorization_bearer':[b'authorization'],
            'credential_assignment':[b'api',b'access',b'auth',b'password',b'passwd',b'secret']}
        for rule,pattern in SECRET_RULES.items():
            if not any(word in lower for word in probes[rule]):continue
            for match in pattern.finditer(data):
                findings.append(dict(path=label,rule=rule,line=data.count(b'\n',0,match.start())+1,
                    match_sha256=hashlib.sha256(match.group()).hexdigest()))
        # Detect common complete environment exports without disclosing any values.
        if all(name in data for name in [b'HOME',b'PATH',b'USER',b'SHELL',b'PWD']) and all(re.search(rb'(?m)["\x27]?'+name+rb'["\x27]?\s*[:=]',data) for name in [b'HOME',b'PATH',b'USER',b'SHELL',b'PWD']):
            findings.append(dict(path=label,rule='possible_complete_environment_dump',line=0,match_sha256=hashlib.sha256(data).hexdigest()))
    for label,p in sorted(paths.items()):
        p=Path(p)
        if p.is_symlink():continue
        data=p.read_bytes();file_hashes[label]=hashlib.sha256(data).hexdigest();scan(data,label)
        try:
            members=[]
            if p.name.endswith(('.tar.gz','.tgz','.tar.bz2','.tar.xz','.tar')):
                with tarfile.open(p,'r:*') as a:
                    total=0
                    for member in a:
                        if not member.isfile():continue
                        total+=member.size
                        if total>1024**3:raise ValueError('expanded archive exceeds audit limit')
                        scan(a.extractfile(member).read(),label+'!/'+member.name);archive_members+=1
            elif p.suffix=='.zip':
                with zipfile.ZipFile(p) as a:
                    if sum(x.file_size for x in a.infolist())>1024**3:raise ValueError('expanded archive exceeds audit limit')
                    for member in a.infolist():
                        if not member.is_dir():scan(a.read(member),label+'!/'+member.filename);archive_members+=1
        except (tarfile.TarError,zipfile.BadZipFile,OSError,ValueError) as e:
            archive_errors.append(dict(path=label,error=type(e).__name__))
    approved={(x['path'],x['rule'],x['line'],x['match_sha256']):x['reason'] for x in CONFIG.get('reviewed_false_positives',[])}
    unresolved=[]
    for finding in findings:
        key=tuple(finding[k] for k in ('path','rule','line','match_sha256'))
        if key in approved:finding['reviewed_false_positive']=approved[key]
        else:unresolved.append(finding)
    return dict(status='PASS' if not unresolved and not archive_errors else 'REVIEW_REQUIRED',scanned_files_and_members=count,
        archive_members=archive_members,findings=findings,unresolved=unresolved,archive_errors=archive_errors,file_sha256=file_hashes,
        method='Credential filename exclusions, byte patterns including logs and metadata, expanded archive member scan; values never written to report',
        limitation='Static pattern checks are not a mathematical proof that arbitrary encoded secrets are absent. Random seeds and ordinary paths retained.')

def links_audit(items,source=None):
    rows=[]
    for key,item in sorted(items.items()):
        p=source/key if source else Path(item['source_path'])
        if not p.is_symlink():continue
        target=os.readlink(p);status='VALID_INTERNAL_RELATIVE'
        try:
            resolved=p.resolve(strict=True)
            root=source if source else ROOTS[Path(key).parts[0]]
            if os.path.isabs(target) or not resolved.is_relative_to(root):status='EXTERNAL_OR_ABSOLUTE'
            elif not source and resolved.relative_to(PROJECTS).as_posix() not in items and not resolved.is_dir():status='TARGET_EXCLUDED'
        except RuntimeError:status='LOOP'
        except OSError:status='BROKEN'
        rows.append(dict(path=key,target=target,status=status))
    return dict(status='PASS' if all(x['status']=='VALID_INTERNAL_RELATIVE' for x in rows) else 'FAIL',links=rows,
        omitted_problem_subtrees=CONFIG['optional_exclusions'],rewritten_links=0)

def process_audit():
    writers=[];inaccessible=0;related=[]
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():continue
        try:
            comm=(proc/'comm').read_text().strip()
            if any(k in comm.lower() for k in ['codex','mirheo','hemocell','mpirun','nvcc','cc1','ninja']):related.append(dict(pid=int(proc.name),comm=comm))
            for fd in (proc/'fd').iterdir():
                try:
                    target=os.readlink(fd)
                    if not any(target.startswith(str(r)+'/') for r in ROOTS.values()):continue
                    flags=re.search(r'^flags:\s+(\d+)',(proc/'fdinfo'/fd.name).read_text(),re.M)
                    if flags and int(flags.group(1),8)&os.O_ACCMODE in (os.O_WRONLY,os.O_RDWR):writers.append(dict(pid=int(proc.name),comm=comm,path=target))
                except (FileNotFoundError,ProcessLookupError):pass
        except (PermissionError,FileNotFoundError,ProcessLookupError):inaccessible+=1
    return dict(related_processes=related,open_source_write_descriptors=writers,inaccessible_or_ended_processes=inaccessible,
        limitation='No process terminated; open-FD check supplemented by per-file before/after stat and SHA-256. Unsaved editor buffers are outside filesystem snapshot.')

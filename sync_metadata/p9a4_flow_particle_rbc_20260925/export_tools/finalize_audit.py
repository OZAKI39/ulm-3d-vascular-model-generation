from pathlib import Path
import os,json,gzip,hashlib,subprocess,shutil,collections,re,datetime
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot']
def sha(p,compressed=False):
 with (gzip.open(p,'rb') if compressed else p.open('rb')) as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
# Preserve final test logs independently of original scientific reports.
shutil.copytree(A/'portable_checks_final',M/'validation',dirs_exist_ok=True)
# Check every original worktree without Git refreshing any original index.
before=json.loads((M/'source_git_before.json').read_text());states=[]
for b in before:
 def git(*args):return subprocess.check_output(['git',*args],cwd=b['root'],env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'}).decode().rstrip('\n')
 s=dict(root=b['root'],head=git('rev-parse','HEAD'),branch=git('branch','--show-current'),index_sha256=sha(Path(b['index_path'])),status=git('status','--porcelain=v1','--untracked-files=all'))
 s['head_unchanged']=s['head']==b['head'];s['branch_unchanged']=s['branch']==b['branch'];s['index_unchanged']=s['index_sha256']==b['index_sha256'];s['status_unchanged']=s['status'].strip()==b['status'].strip();states.append(s)
dump(M/'source_git_after.json',states)
# Verify copied P9 scientific inputs, outputs and code against original bytes.
local=json.loads((M/'local_added_inventory.json').read_text());packed={r['source']:r for r in json.loads((M/'packed_artifacts.json').read_text())};p9=[]
for r in local:
 if r['group']!='P9A4':continue
 source=Path(r['source']);p=packed.get(r['source']);target=D/(p['packed_path'] if p else r['target']);expected=p['original_sha256'] if p else r['sha256']
 p9.append(dict(source=str(source),target=str(target.relative_to(D)),expected_sha256=expected,original_sha256=sha(source),published_original_sha256=sha(target,p is not None)))
assert len(p9)==321 and all(r['original_sha256']==r['expected_sha256']==r['published_original_sha256'] for r in p9)
dump(M/'p9a4_byte_validation.json',dict(files=len(p9),all_original_and_published_bytes_identical=True,files_detail=p9))
# Record whether independently edited vascular files changed again after snapshot.
v=json.loads((M/'vascular_source_refresh.json').read_text());changes=[]
for r in v['files']:
 p=Path(r['source']);actual=sha(p) if p.is_file() else None
 if actual!=r['sha256']:changes.append(dict(source=r['source'],copied_sha256=r['sha256'],current_sha256=actual))
dump(M/'concurrent_source_changes.json',dict(checked_utc=now(),vascular_files_checked=len(v['files']),vascular_changes_since_copy=changes,policy='Source trees are read-only; a captured version is not overwritten by concurrent development.'))
# Re-resolve SHA references after navigation-only documentation edits.
x=json.loads(gzip.decompress((M/'server_delta_inventory.json.gz').read_bytes()));cache={};fixed=[];invalid=[];exported=0
for r in x['files']:
 if r['status']=='EXPORTED':
  p=D/r['target'];actual=sha(p)
  if r['representation']=='IDENTICAL_BYTES':
   if actual!=r['sha256']:invalid.append(dict(source=r['source'],target=r['target'],reason='EXPORTED_SHA_MISMATCH'))
  elif r['representation']=='POINTER_EVIDENCE_JSON':
   if not isinstance(json.loads(p.read_text()),dict):raise ValueError('Invalid pointer receipt')
  exported+=1
 elif r['status']=='REPRESENTED_BY_EXISTING_ARTIFACT':
  rep=r['representation'];key=(rep['path'],rep['representation']);p=D/rep['path']
  if key not in cache:cache[key]=sha(p,rep['representation']=='GZIP_OF_IDENTICAL_BYTES') if p.is_file() else None
  if cache[key]!=r['sha256']:
   saved=M/'original_docs'/rep['path']
   if saved.is_file() and sha(saved)==r['sha256']:
    fixed.append(dict(source=r['source'],old_target=rep['path'],new_target=str(saved.relative_to(D))))
    rep['path']=str(saved.relative_to(D));rep['representation']='IDENTICAL_BYTES'
   else:invalid.append(dict(source=r['source'],target=rep['path'],reason='CONTENT_REFERENCE_SHA_MISMATCH'))
assert not invalid,invalid[:5]
if fixed:
 with (M/'server_delta_inventory.json.gz').open('wb') as raw:
  with gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0,compresslevel=6) as f:f.write((json.dumps(x,indent=2)+'\n').encode())
dump(M/'server_reference_validation.json',dict(checked_utc=now(),exported_files_verified=exported,reference_records_verified=sum(r['status']=='REPRESENTED_BY_EXISTING_ARTIFACT' for r in x['files']),unique_targets=len(cache),navigation_reference_remaps=fixed,invalid=invalid))
s=json.loads((M/'server_delta_summary.json').read_text());s['inventory_compressed_bytes']=(M/'server_delta_inventory.json.gz').stat().st_size;dump(M/'server_delta_summary.json',s)
# Scan newly published/updated bytes, without exposing any matched value.
old=set(subprocess.check_output(['git','ls-files','-z'],cwd=D).decode().split('\0'));changed=set(subprocess.check_output(['git','diff','--name-only','-z'],cwd=D).decode().split('\0'))
allfiles=sorted(p for p in D.rglob('*') if p.is_file() and p.name!='.git');newfiles=[p for p in allfiles if str(p.relative_to(D)) not in old or str(p.relative_to(D)) in changed]
patterns={
 'private_key':rb'-----BEGIN (?:RSA |EC |DSA |OPENSSH |ENCRYPTED )?PRIVATE KEY-----',
 'github_token':rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b',
 'aws_access_key':rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b',
 'openai_project_key':rb'\bsk-(?:proj-|svcacct-)[A-Za-z0-9_-]{30,}\b',
 'slack_token':rb'\bxox[baprs]-[A-Za-z0-9-]{20,}\b',
 'credentialed_url':rb'https?://[^\s/@:]{2,}:[^\s/@]{8,}@(?!localhost|127\.0\.0\.1)',
}
findings=[];pointers=[];symlinks=[];suspect_names=[]
for p in newfiles:
 rel=str(p.relative_to(D))
 if p.is_symlink():symlinks.append(rel)
 if p.name=='.env' or (p.suffix in {'.key','.p12','.pfx'} or p.name in {'id_rsa','id_ed25519','credentials.json'}):suspect_names.append(rel)
 raw=p.read_bytes()
 if raw.startswith(b'version https://git-lfs.github.com/spec/v1\n'):pointers.append(rel)
 if p.suffix=='.gz':
  try:raw=gzip.decompress(raw)
  except (OSError,EOFError):continue
 for name,pattern in patterns.items():
  if re.search(pattern,raw):findings.append(dict(path=rel,pattern=name,sha256=sha(p)))
result=dict(checked_utc=now(),new_or_modified_files_scanned=len(newfiles),scope='New/modified snapshot bytes and decompressed gzip content; obvious credential signatures and private-key headers; not a proof that every possible secret format is absent.',findings=findings,suspect_filenames=suspect_names,lfs_pointers=pointers,symlinks=symlinks,total_files=len(allfiles),total_bytes=sum(p.stat().st_size for p in allfiles),largest_files=[dict(path=str(p.relative_to(D)),bytes=p.stat().st_size) for p in sorted(allfiles,key=lambda p:p.stat().st_size,reverse=True)[:10]])
dump(M/'preflight_audit.json',result)
summary=dict(p9_scientific_byte_checks=len(p9),source_git=[{k:s[k] for k in ['root','head_unchanged','branch_unchanged','index_unchanged','status_unchanged']} for s in states],concurrent_vascular_changes=len(changes),server_exported_files_verified=exported,server_reference_remaps=len(fixed),scan_findings=findings,suspect_filenames=suspect_names,lfs_pointers=pointers)
print(json.dumps(summary,indent=2));assert not findings and not suspect_names and not pointers and not symlinks

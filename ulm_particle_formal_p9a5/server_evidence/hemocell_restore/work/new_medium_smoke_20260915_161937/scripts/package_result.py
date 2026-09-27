from pathlib import Path
import tarfile,hashlib,subprocess,json,time
W=Path('/workspace/hemocell_restore/work/new_medium_smoke_20260915_161937');F=Path('/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937');D=W/'reports/transport';D.mkdir(parents=True,exist_ok=False);out=D/'new_medium_smoke_result.tar.zst'
expected={name:digest for digest,name in (line.split('  ',1) for line in (F/'SHA256SUMS').read_text().splitlines())};expected['SHA256SUMS']=hashlib.sha256((F/'SHA256SUMS').read_bytes()).hexdigest()
seen={};links=0;regular=0;t=time.monotonic()
proc=subprocess.Popen(['zstd','-q','-T2','-9','-o',str(out)],stdin=subprocess.PIPE)
with tarfile.open(fileobj=proc.stdin,mode='w|',format=tarfile.PAX_FORMAT) as tar:
 for p in sorted(F.rglob('*')):
  rel=str(p.relative_to(F));info=tar.gettarinfo(str(p),arcname=rel)
  if p.is_dir():tar.addfile(info);continue
  assert p.is_file() and rel in expected
  with p.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
  assert actual==expected[rel],rel
  key=(actual,info.mode,info.uid,info.gid)
  if key in seen:
   info.type=tarfile.LNKTYPE;info.linkname=seen[key];info.size=0;tar.addfile(info);links+=1
  else:
   with p.open('rb') as f:tar.addfile(info,f)
   seen[key]=rel;regular+=1
proc.stdin.close();assert proc.wait()==0
with out.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
r=dict(status='PASS',path=str(out),bytes=out.stat().st_size,sha256=digest,root_manifest_sha256=expected['SHA256SUMS'],file_count=len(expected),unique_regular_files=regular,lossless_hardlink_entries=links,source_files_changed=False,seconds=time.monotonic()-t,scope='Transport only: exact file bytes, names and modes restored; duplicate identical regular files represented as hardlinks')
(D/'TRANSPORT_PACKAGE.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

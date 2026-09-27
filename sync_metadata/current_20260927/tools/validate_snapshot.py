"""Verify copied bytes, protected inputs, internal links and publication inventory."""
from pathlib import Path
import json,csv,hashlib,re,collections
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260927'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 failures=[];checked=0
 for row in csv.DictReader((M/'local_file_inventory.csv').open()):
  src=Path(row['source']);dst=D/row['destination']
  if not src.is_file() or not dst.is_file() or sha(src)!=row['sha256'] or sha(dst)!=row['sha256']:failures.append(dict(type='local_bytes',source=str(src),snapshot=str(dst)))
  checked+=1
 server=json.loads((M/'server_file_map.json').read_text());hash_cache={}
 pointers=json.loads((M/'upstream_lfs_pointers.json').read_text())
 pointer_hashes={r['sha256']:hashlib.sha256(r['content'].encode()).hexdigest() for r in pointers}
 for r in pointers:
  if sha(Path(r['source']))!=r['sha256'] or pointer_hashes[r['sha256']]!=r['sha256']:failures.append(dict(type='archived_lfs_pointer',source=r['source']))
 for r in server:
  if r['disposition']=='upstream_lfs_pointer_metadata_only':
   if pointer_hashes.get(r['sha256'])!=r['sha256']:failures.append(dict(type='archived_server_lfs_pointer',source=r['source']))
   continue
  p=D/r['snapshot']
  if str(p) not in hash_cache:hash_cache[str(p)]=sha(p) if p.is_file() else None
  if hash_cache[str(p)]!=r['sha256']:failures.append(dict(type='server_bytes',source=r['source'],snapshot=r['snapshot']))
 links=[];files=[];pointers=[];security=[]
 patterns={
  'private_key':rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----',
  'github_token':rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b',
  'aws_access_key':rb'\b(?:AKIA|ASIA)[0-9A-Z]{16}\b',
  'slack_token':rb'\bxox[baprs]-[0-9A-Za-z-]{20,}\b',
  'openai_key':rb'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{40,}\b',
  'credential_in_url':rb'https?://[^\s/:@]{2,}:[^\s/@]{6,}@',
 }
 triggers={'private_key':[b'-----BEGIN '],'github_token':[b'ghp_',b'gho_',b'ghu_',b'ghs_',b'ghr_',b'github_pat_'],'aws_access_key':[b'AKIA',b'ASIA'],'slack_token':[b'xox'],'openai_key':[b'sk-'],'credential_in_url':[b'http://',b'https://']}
 for p in sorted(D.rglob('*')):
  rel=str(p.relative_to(D))
  if rel=='.git':continue
  if p.is_symlink():
   ok=p.exists() and p.resolve().is_relative_to(D)
   links.append(dict(path=rel,target=str(p.readlink()),valid=ok))
   if not ok:failures.append(dict(type='symlink',snapshot=rel))
   continue
  if not p.is_file():continue
  payload=p.read_bytes();n=len(payload)
  files.append(dict(path=rel,bytes=n,sha256=hashlib.sha256(payload).hexdigest()))
  if n>=50*1024**2:failures.append(dict(type='oversize',snapshot=rel,bytes=n))
  if payload.startswith(b'version https://git-lfs.github.com/spec/v1\n'):
   pointers.append(dict(path=rel,content=payload.decode()))
  # Only record path/line/category, never the potential credential value.
  if rel.startswith('sync_metadata/current_20260927/validation/'):continue
  for kind,pat in patterns.items():
   if not any(t in payload for t in triggers[kind]):continue
   for match in re.finditer(pat,payload):
    security.append(dict(path=rel,line=payload[:match.start()].count(b'\n')+1,kind=kind))
 result=dict(local_file_byte_checks=checked,remote_mapped_file_checks=len(server),files=len(files),bytes=sum(x['bytes'] for x in files),symlinks=len(links),failures=failures,lfs_pointer_files=len(pointers),credential_pattern_candidates=len(security))
 (M/'validation/snapshot_validation.json').write_text(json.dumps(result,indent=2)+'\n')
 (M/'validation/lfs_pointer_inventory.json').write_text(json.dumps(pointers,indent=2)+'\n')
 (M/'validation/security_pattern_candidates.json').write_text(json.dumps(security,indent=2)+'\n')
 (M/'validation/symlink_validation.json').write_text(json.dumps(links,indent=2)+'\n')
 with (M/'snapshot_manifest.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);w.writeheader()
  for r in files:
   if r['path'].startswith('sync_metadata/current_20260927/validation/') or r['path']=='sync_metadata/current_20260927/snapshot_manifest.csv':continue
   w.writerow(r)
 print(json.dumps(result,indent=2))
 if failures:raise SystemExit(1)
if __name__=='__main__':main()

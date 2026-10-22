from pathlib import Path
import json,hashlib,os,time
root=Path('/workspace/vascular_printing');out=root/'migration/20261009';start=time.monotonic()
rows=json.loads((out/'project_files.json').read_text());bad=[];total=0
for name,row in rows.items():
 p=root/name
 if row['type']=='symlink':
  ok=p.is_symlink() and os.readlink(p)==row['target'] and p.exists()
 else:
  h=hashlib.sha256()
  if p.is_file():
   with p.open('rb') as f:
    for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
  ok=p.is_file() and h.hexdigest()==row['sha256']
 if not ok:bad.append(name)
 total+=row['bytes']
external=json.loads((out/'external_files.json').read_text())
for source,row in external.items():
 p=out/'external'/row['target']
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:bad.append(str(p))
r=dict(status='VERIFIED' if not bad else 'FAILED',project_entries=len(rows),external_entries=len(external),source_bytes=total,mismatches=bad,elapsed_seconds=time.monotonic()-start)
(out/'transfer_verification.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
raise SystemExit(bool(bad))

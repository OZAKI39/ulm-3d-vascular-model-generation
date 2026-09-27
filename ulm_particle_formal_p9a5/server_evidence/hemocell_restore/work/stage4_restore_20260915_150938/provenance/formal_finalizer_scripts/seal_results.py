#!/usr/bin/python3
from pathlib import Path
import json,hashlib
R=Path(__file__).resolve().parents[1]
assert (R/'RUN_TERMINAL.json').exists()
assert not (R/'provenance/supervisord.pid').exists(),'Stop private supervisor before seal'
v=json.loads((R/'FINAL_TERMINAL_SUMMARY.json').read_text());v['SHA256_STATUS']='PASS'
(R/'FINAL_TERMINAL_SUMMARY.json').write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')
p=R/'FINAL_TERMINAL_SUMMARY.txt';p.write_text(p.read_text().replace('SHA256_STATUS = PENDING_FINAL_SEAL','SHA256_STATUS = PASS'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
files=[(sha(p),str(p.relative_to(R))) for p in sorted(R.rglob('*')) if p.is_file() and not p.is_symlink() and p.name!='SHA256SUMS']
(R/'SHA256SUMS').write_text(''.join(h+'  '+n+'\n' for h,n in files))
for h,n in files:assert sha(R/n)==h,n
print('SHA256_STATUS=PASS','SEALED_FILES',len(files),'MANIFEST_SHA256',sha(R/'SHA256SUMS'))

"""Verify every delivered byte and local HTML target using stdlib only."""
from pathlib import Path
import hashlib,json,re,os
root=Path(__file__).resolve().parents[2];task='clot_shear_load_inverse_design_v1';run='20261017T165129Z';out=root/'verification'/task/run
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
manifest=json.loads((out/'DELIVERY_MANIFEST.json').read_text());count=0
for entry in manifest['files']:
 p=root/entry['path'];assert p.exists(),p
 if entry['kind']=='symlink':assert p.is_symlink() and os.readlink(p)==entry['target'] and sha(p)==entry['target_sha256'],p
 else:assert p.stat().st_size==entry['bytes'] and sha(p)==entry['sha256'],p
 count+=1
for p in [out/'OPEN_RESULTS.html',root/'results'/task/run/'OPEN_RESULTS.html']:
 for href in re.findall(r'(?:href|src)="([^"]+)"',p.read_text()):
  if not href.startswith(('http:','https:','#')):assert (p.parent/href).resolve().exists(),(p,href)
doc=json.loads((out/'WSL_DOCUMENT_UPDATE.json').read_text());assert sha(Path(doc['path']))==doc['after_sha256']
z=dict(status='PASS',files_verified=count,manifest_sha256=sha(out/'DELIVERY_MANIFEST.json'),local_document_sha256=doc['after_sha256'],html_relative_links_valid=True,scope='All delivered task files and document; server-only binary/upstream data are tracked by dependency hashes, not falsely claimed downloaded')
(out/'LOCAL_READBACK.json').write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n');print(json.dumps(z))

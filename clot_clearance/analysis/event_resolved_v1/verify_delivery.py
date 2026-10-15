"""Verify sealed payload on either server or WSL; never run dynamics."""
from pathlib import Path
import hashlib,json,sys
r=Path(sys.argv[1]);task='lps_native_softening_event_resolved_integration_v1';run='20261013T235945Z';p=r/'verification'/task/run
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
manifest=p/'PAYLOAD_MANIFEST.json';assert sha(manifest)==(p/'PAYLOAD_MANIFEST.sha256').read_text().strip()
z=json.loads(manifest.read_text());errors=[]
for name,m in z['files'].items():
 q=r/name
 if not q.is_file() or q.stat().st_size!=m['bytes'] or sha(q)!=m['sha256']:errors.append(name)
assert not errors,errors
doc=json.loads((p/'WSL_DOCUMENT_UPDATE.json').read_text())
if str(r).startswith('/home/'):
 assert sha(Path(doc['actual_path']))==doc['readback_sha256']
 assert (p/'wsl_document/clot_rupture.after.md').read_bytes()==Path(doc['actual_path']).read_bytes()
receipt=dict(status='ALL_SEALED_FILES_SHA256_VERIFIED',root=str(r),file_count=z['file_count'],total_bytes=z['total_bytes'],manifest_sha256=sha(manifest),actual_WSL_document_checked=str(r).startswith('/home/'),full700_status='FAILED_ATOM_LOSS',local_event='PASS')
(p.parent/('LOCAL_SYNC_RECEIPT_'+run+'.json' if str(r).startswith('/home/') else 'SERVER_VERIFICATION_RECEIPT_'+run+'.json')).write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))

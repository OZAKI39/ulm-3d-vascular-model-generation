"""Verify the selected GitHub snapshot with the Python standard library only."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
manifest=json.loads((Path(__file__).resolve().parent/'SYNC_FILE_MANIFEST.json').read_text())
errors=[]
for row in manifest['files']:
 p=ROOT/row['path']
 actual=hashlib.file_digest(p.open('rb'),'sha256').hexdigest() if p.is_file() else None
 if actual!=row['sha256']:errors.append(row['path'])
assert not errors,errors
print(json.dumps(dict(status='PASS',verified_files=len(manifest['files']),bytes=sum(r['size'] for r in manifest['files']))))

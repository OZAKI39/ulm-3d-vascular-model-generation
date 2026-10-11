"""Verify this sync's explicit payload; historical full-run manifests are separate."""
from pathlib import Path
import hashlib,json

def main():
 here=Path(__file__).resolve().parent;root=here.parent
 manifest=json.loads((here/'SNAPSHOT_MANIFEST.json').read_text());failures=[]
 for row in manifest['files']:
  path=root/row['path']
  if not path.resolve().is_relative_to(root.resolve()) or not path.is_file():
   failures.append(row['path']+': missing/unsafe');continue
  with path.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
  if actual!=row['sha256'] or path.stat().st_size!=row['bytes']:failures.append(row['path']+': checksum/size')
 if failures:raise RuntimeError('\n'.join(failures))
 print(json.dumps(dict(status='PASS',files=len(manifest['files']),bytes=sum(r['bytes'] for r in manifest['files']),historical_files_not_reinterpreted=True)))

if __name__=='__main__':main()

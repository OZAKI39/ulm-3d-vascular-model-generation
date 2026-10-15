"""Check this incremental GitHub snapshot, not the omitted full-server payload."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];G=R/'github_sync_20261015'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def main():
 m=json.loads((G/'SNAPSHOT_MANIFEST.json').read_text())
 for row in m['files']:
  p=R/row['path'];assert p.resolve().is_relative_to(R.resolve())
  assert p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],str(p)
 print(json.dumps(dict(status='PASS',files=len(m['files']),bytes=sum(x['bytes'] for x in m['files']),base_commit=m['base_commit']),indent=2))
if __name__=='__main__':main()

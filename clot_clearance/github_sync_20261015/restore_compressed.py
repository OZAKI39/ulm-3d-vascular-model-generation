"""Verify or restore the lossless text exports without replacing differing files."""
from pathlib import Path
import argparse,gzip,hashlib,json
R=Path(__file__).resolve().parents[1]
def digest(b):return hashlib.sha256(b).hexdigest()
def entries():return json.loads((R/'github_sync_20261015/SOURCE_FILES.json').read_text())['files']
def unpack(row):
 p=R/row['path'];b=p.read_bytes()
 assert digest(b)==row['sha256'],str(p)
 raw=gzip.decompress(b)
 assert digest(raw)==row['uncompressed_export_sha256'],str(p)
 if not row['credential_redactions']:assert digest(raw)==row['source_sha256'],str(p)
 return raw
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify-only',action='store_true');args=parser.parse_args();count=0
 for row in entries():
  if row['compression']!='gzip':continue
  raw=unpack(row);p=R/row['source_path'];assert p.resolve().is_relative_to(R.resolve())
  if not args.verify_only:
   if p.exists():assert p.read_bytes()==raw,'Refusing to overwrite differing file: '+str(p)
   else:
    with p.open('xb') as f:f.write(raw)
  count+=1
 print(json.dumps(dict(status='PASS',compressed_files=count,mode='verify_only' if args.verify_only else 'restore'),indent=2))
if __name__=='__main__':main()

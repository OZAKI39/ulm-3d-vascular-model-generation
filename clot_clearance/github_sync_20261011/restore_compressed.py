"""Verify or restore exact scientific CSV/log bytes from the reviewed export."""
from pathlib import Path
import argparse,gzip,hashlib,json,os

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--verify-only',action='store_true')
 args=parser.parse_args();here=Path(__file__).resolve().parent;root=here.parent
 rows=json.loads((here/'SOURCE_FILES.json').read_text())['files'];count=0
 for row in rows:
  if row['compression']!='gzip_lossless':continue
  packed=root/row['path'];raw=gzip.decompress(packed.read_bytes())
  if hashlib.sha256(raw).hexdigest()!=row['source_sha256'] or len(raw)!=row['source_bytes']:
   raise RuntimeError('Source-byte verification failed: '+row['path'])
  target=root/row['source_path']
  if not target.resolve().is_relative_to(root.resolve()):raise RuntimeError('Unsafe output path')
  if not args.verify_only:
   if target.exists():
    if hashlib.sha256(target.read_bytes()).hexdigest()!=row['source_sha256']:
     raise RuntimeError('Refusing to replace different existing data: '+str(target))
   else:
    target.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive create: another process cannot be overwritten after the existence check.
    with target.open('xb') as f:f.write(raw)
  count+=1
 print(json.dumps(dict(status='PASS',verified_compressed_files=count,restored=not args.verify_only)))

if __name__=='__main__':main()

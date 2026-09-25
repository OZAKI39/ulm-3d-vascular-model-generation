"""Restore exact P9-A.4 scientific bytes from gzip; never replace existing files."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,shutil,tempfile
ROOT=Path(__file__).resolve().parents[2]
META=Path(__file__).resolve().parent

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def restore(destination,verify_only=False):
 destination=Path(destination).resolve();rows=json.loads((META/'packed_artifacts.json').read_text());result=[]
 for row in rows:
  packed=(ROOT/row['packed_path']).resolve();target=(destination/row['original_path']).resolve()
  if not packed.is_relative_to(ROOT) or not target.is_relative_to(destination):raise ValueError('Artifact escapes root')
  if sha(packed)!=row['packed_sha256']:raise ValueError('Packed SHA mismatch: '+str(packed))
  with gzip.open(packed,'rb') as f:
   if hashlib.file_digest(f,'sha256').hexdigest()!=row['original_sha256']:raise ValueError('Decompressed SHA mismatch')
  status='VERIFIED_COMPRESSED_AND_ORIGINAL_BYTES'
  if not verify_only:
   if target.exists():
    if sha(target)!=row['original_sha256']:raise FileExistsError('Refusing to overwrite different bytes: '+str(target))
    status='EXISTING_IDENTICAL'
   else:
    target.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix='.p9a4-restore-',dir=target.parent,delete=False) as tmp:
     name=Path(tmp.name)
     with gzip.open(packed,'rb') as f:shutil.copyfileobj(f,tmp,1024**2)
     tmp.flush();os.fsync(tmp.fileno())
    try:
     if sha(name)!=row['original_sha256']:raise ValueError('Restored bytes mismatch')
     os.link(name,target)  # atomic exclusive publication; fails if another writer won
    finally:name.unlink()
    status='RESTORED'
  result.append(dict(path=row['original_path'],status=status,sha256=row['original_sha256']))
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-only',action='store_true');p.add_argument('--destination',type=Path,default=ROOT);a=p.parse_args()
 print(json.dumps(restore(a.destination,a.verify_only),indent=2))

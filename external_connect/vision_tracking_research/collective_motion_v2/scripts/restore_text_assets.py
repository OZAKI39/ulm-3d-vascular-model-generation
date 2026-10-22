"""Restore compressed numeric records in a fresh Git checkout, never overwrite."""
from pathlib import Path
import gzip,hashlib
ROOT=Path(__file__).resolve().parents[1]
for file in (ROOT/'results').rglob('*.gz'):
 target=file.with_suffix('');data=gzip.decompress(file.read_bytes())
 if target.exists():
  if hashlib.sha256(target.read_bytes()).digest()!=hashlib.sha256(data).digest():raise RuntimeError('DIFFERENT_EXISTING_FILE:'+str(target))
 else:target.write_bytes(data)
print('Text assets available')

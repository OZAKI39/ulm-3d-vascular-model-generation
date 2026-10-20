"""Restore losslessly compressed input evidence in a compact Git checkout."""
from pathlib import Path
import json,gzip,hashlib,os
R=Path(__file__).resolve().parents[3];O=R/'verification/clot_irregular_lps_discretization_v4/20261019T221000Z'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 spec=json.loads((O/'GIT_COMPRESSED_EVIDENCE.json').read_text())
 for q in spec['records']:
  compressed=O/q['compressed'];target=O/q['source'];assert target.resolve().is_relative_to(O.resolve());raw=compressed.read_bytes();assert sha(raw)==q['compressed_sha256'];data=gzip.decompress(raw);assert len(data)==q['source_bytes'] and sha(data)==q['source_sha256']
  if target.exists():assert target.read_bytes()==data,f'Refusing to overwrite different evidence: {target}'
  else:
   target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_name(target.name+'.restore.tmp');tmp.write_bytes(data);os.replace(tmp,target)
  print('VERIFIED',q['source'])
if __name__=='__main__':main()

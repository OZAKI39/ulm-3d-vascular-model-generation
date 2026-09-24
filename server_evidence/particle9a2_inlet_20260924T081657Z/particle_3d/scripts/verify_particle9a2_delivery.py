#!/usr/bin/env python3
"""Hash all deliverables and independently verify a synchronized server copy."""
import argparse,json,hashlib,subprocess
from pathlib import Path
from particle_3d.injection_method_c import canonical_bytes

ROOT=Path(__file__).resolve().parents[2];R=ROOT/'particle_3d/reports/particle9a2_inlet_sampling';D=R/'data'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--verify-remote',action='store_true');a=p.parse_args()
    if not a.verify_remote:
        files=list(R.rglob('*'))+list((ROOT/'particle_3d/outputs/particle9a2_2mmps').rglob('*'))
        source=json.loads((D/'new_source_manifest.json').read_text());files.extend(ROOT/k for k in source)
        files.extend(ROOT/'particle_3d/contracts'/n for n in ['SONOVUE_DIAMETER_MAX4UM_V1.json','PARTICLE9A2_PRODUCTION_V1.json'])
        exempt={'delivery_manifest.json','server_delivery_verification.json'}
        hashes={str(f.relative_to(ROOT)):sha(f) for f in sorted(set(files)) if f.is_file() and f.name not in exempt and '__pycache__' not in f.parts}
        manifest=dict(schema='P9A2_DELIVERY_SHA256_V1',files=hashes,self_reference_exemptions=sorted(exempt))
        (D/'delivery_manifest.json').write_bytes(canonical_bytes(manifest));print('manifest',len(hashes));return
    remote=json.loads((D/'remote.json').read_text())['root'];manifest=json.loads((D/'delivery_manifest.json').read_text())
    for rel,h in manifest['files'].items():assert sha(ROOT/rel)==h,rel
    program='''from pathlib import Path
import json,hashlib
root=Path(ROOT_LITERAL)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
m=json.loads((root/'particle_3d/reports/particle9a2_inlet_sampling/data/delivery_manifest.json').read_text())
bad=[k for k,v in m['files'].items() if not (root/k).is_file() or sha(root/k)!=v]
print(json.dumps(dict(status='MATCH' if not bad else 'MISMATCH',file_count=len(m['files']),mismatches=bad,server_root=str(root),manifest_sha256=sha(root/'particle_3d/reports/particle9a2_inlet_sampling/data/delivery_manifest.json'))))
'''.replace('ROOT_LITERAL',repr(remote))
    cmd=['ssh','-p','4159','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@50.115.148.16','python3 -']
    proc=subprocess.run(cmd,input=program,text=True,capture_output=True,check=True);result=json.loads(proc.stdout)
    assert result['status']=='MATCH',result
    (D/'server_delivery_verification.json').write_bytes(canonical_bytes(result));print(json.dumps(result))
if __name__=='__main__':main()

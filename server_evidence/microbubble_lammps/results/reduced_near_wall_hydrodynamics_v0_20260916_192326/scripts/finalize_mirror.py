"""Mirror only this new stage and verify every artifact by SHA256."""
import csv
import hashlib
import io
import json
import subprocess
import tarfile
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
paths=json.loads((ROOT/'provenance/STAGE_PATHS.json').read_text())
ssh=['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=yes',
     '-o','UpdateHostKeys=no','-o','ControlMaster=no','-o','ControlPath=none','vast4090']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def files():return sorted(p for p in ROOT.rglob('*') if p.is_file())
def upload(items):
    buffer=io.BytesIO()
    with tarfile.open(fileobj=buffer,mode='w:gz') as tar:
        for p in items:tar.add(p,arcname=str(p.relative_to(ROOT)))
    p=subprocess.run(ssh+['tar -xzf - -C '+paths['remote_stage']],input=buffer.getvalue(),capture_output=True,timeout=60)
    assert p.returncode==0,p.stderr
def verify(expected):
    code='''from pathlib import Path
import hashlib,json
r=Path(__ROOT__);expected=__EXPECTED__
actual={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in r.rglob('*') if p.is_file()}
assert actual==expected,'REMOTE_LOCAL_PAYLOAD_MISMATCH'
print(json.dumps({'match':True,'checked_files':len(actual),'sha256':actual}))
'''.replace('__ROOT__',repr(paths['remote_stage'])).replace('__EXPECTED__',repr(expected))
    p=subprocess.run(ssh+['python3 -B -'],input=code,text=True,capture_output=True,timeout=45)
    assert p.returncode==0,p.stderr
    return json.loads(p.stdout)
def manifests():
    source=[p for p in files() if p.suffix in ['.py','.cpp','.hpp'] or p.name=='CMakeLists.txt']
    with (ROOT/'provenance/SOURCE_SHA256.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['relative_path','sha256','size'])
        w.writerows((str(p.relative_to(ROOT)),sha(p),p.stat().st_size) for p in source)
    payload=[p for p in files() if p.name not in ['PAYLOAD_SHA256.tsv','REMOTE_LOCAL_SHA_CHECK.json']]
    with (ROOT/'provenance/PAYLOAD_SHA256.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['relative_path','sha256','size'])
        w.writerows((str(p.relative_to(ROOT)),sha(p),p.stat().st_size) for p in payload)
    return {str(p.relative_to(ROOT)):sha(p) for p in files()}

assert json.loads((ROOT/'provenance/IMMUTABILITY_CHECK.json').read_text())['previous_remote_stage_unchanged']
assert json.loads((ROOT/'validation/PHASE1_VALIDATION.json').read_text())['ALGEBRA_PASS']
first=manifests();upload(files());verify(first)
p=ROOT/'CURRENT_STATE.json';state=json.loads(p.read_text());state['mirror_status']='VERIFIED_SHA256'
state['source_remote_local_sha_match']=True;p.write_text(json.dumps(state,indent=2)+'\n')
final=manifests();upload(files());checked=verify(final)
receipt={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'local_stage':str(ROOT),'remote_stage':paths['remote_stage'],
         'payload_match':True,'source_match':True,'payload_files_checked_excluding_receipt':len(final),
         'final_files_including_receipt':len(final)+1,'payload_bytes_excluding_receipt':sum(p.stat().st_size for p in files()),
         'sha256':checked['sha256'],'receipt_self_hash_note':'Finishing routine separately verifies the uploaded receipt hash; receipt is excluded from its own recursive manifest',
         'build_objects_mirrored':False}
p=ROOT/'provenance/REMOTE_LOCAL_SHA_CHECK.json';p.write_text(json.dumps(receipt,indent=2)+'\n');upload([p])
final[str(p.relative_to(ROOT))]=sha(p);last=verify(final)
print(json.dumps({'FINAL_MIRROR_VERIFIED':True,'files':last['checked_files'],'receipt_sha256':sha(p),'source_remote_local_match':True}))
print((ROOT/'validation/FINAL_STATUS.txt').read_text())

#!/usr/bin/env python3
"""Collect immutable remote stage outputs and verify every transported byte."""
from pathlib import Path
import argparse,hashlib,json,shlex,subprocess
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2]
REMOTE='/workspace/particle8_2a_20260922/results'
SSH=['ssh','-p','4159','-i',str(Path.home()/'.ssh/vast_step3b_ed25519'),'-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@50.115.148.16']


def main():
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');a=p.parse_args()
    out=ROOT/'particle_3d/outputs/particle8_2a';out.mkdir(parents=True,exist_ok=True)
    report=ROOT/'particle_3d/reports/particle8_2a'
    if not a.partial:
        code='''from pathlib import Path
import hashlib,json,socket,time
r=Path('/workspace/particle8_2a_20260922/results')
assert (r/'COMPUTE_AND_MEDIA_COMPLETE').exists()
assert (r/'FINAL_QUALITY_COMPLETE').exists()
assert (r/'review/REMOTE_IMMUTABILITY_VERIFICATION.json').exists()
(r/'STOP_MONITOR').touch()
time.sleep(11)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
files={str(p.relative_to(r)):digest(p) for p in sorted(r.rglob('*')) if p.is_file() and p.name!='REMOTE_ARTIFACT_SHA256.json'}
(r/'REMOTE_ARTIFACT_SHA256.json').write_text(json.dumps(dict(hostname=socket.gethostname(),sha256=files),indent=2)+'\\n')
print('REMOTE_ARTIFACTS',len(files))
'''
        subprocess.run(SSH+['python3 -c '+shlex.quote(code)],check=True)
    subprocess.run(['rsync','-a','--partial','--delay-updates','--info=stats1','-e',shlex.join(SSH[:-1]),SSH[-1]+':'+REMOTE+'/',str(out)+'/'],check=True)
    if a.partial:return
    inventory=json.loads((out/'REMOTE_ARTIFACT_SHA256.json').read_text());mismatches=[]
    for name,expected in inventory['sha256'].items():
        h=hashlib.sha256()
        with (out/name).open('rb') as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
        if h.hexdigest()!=expected:mismatches.append(name)
    receipt=dict(all_pass=not mismatches,hostname=inventory['hostname'],file_count=len(inventory['sha256']),
        mismatches=mismatches,verified_utc=datetime.now(timezone.utc).isoformat(),
        remote_manifest_sha256=hashlib.sha256((out/'REMOTE_ARTIFACT_SHA256.json').read_bytes()).hexdigest(),
        remote_root=REMOTE,local_root=str(out),method='SHA256_EVERY_TRANSFERRED_REMOTE_ARTIFACT')
    (report/'REMOTE_LOCAL_CHECKSUM_VERIFICATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    assert not mismatches
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()

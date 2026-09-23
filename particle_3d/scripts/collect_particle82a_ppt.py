#!/usr/bin/env python3
"""Collect only the user-revised 500-track presentation and verify its provenance."""
from pathlib import Path
import hashlib,json,shlex,subprocess

ROOT=Path(__file__).resolve().parents[2]
REMOTE='/workspace/particle8_2a_20260922'
COMMITS=['83c0cc3ccb5cbc88865ba4e3dcbeb92438f299bd','4a38c85a698c21f4188ed345a7a0b34e947193b6']
SSH=['ssh','-p','4159','-i',str(Path.home()/'.ssh/vast_step3b_ed25519'),'-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@50.115.148.16']


def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    out=ROOT/'particle_3d/outputs/particle8_2a_ppt';out.mkdir(exist_ok=True)
    report=ROOT/'particle_3d/reports/particle8_2a_ppt';report.mkdir(exist_ok=True)
    code=r'''from pathlib import Path
import hashlib,json,socket
base=Path(%r);out=base/'results/ppt'
assert (out/'PPT_COMPLETE').exists(), 'Presentation render/audit/tests not complete'
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
records=[]
for commit in %r:
 p=base/commit;manifest=json.loads((p/'provenance.json').read_text());counts={}
 for group in ['source_sha256','frozen_input_sha256','external_scientific_input_sha256']:
  for file,expected in manifest[group].items():assert digest(p/'repo'/file)==expected,file
  counts[group]=len(manifest[group])
 records.append(dict(source_commit=commit,counts=counts,provenance_sha256=digest(p/'provenance.json')))
check=dict(all_pass=True,hostname=socket.gethostname(),snapshots=records)
(out/'PPT_REMOTE_IMMUTABILITY.json').write_text(json.dumps(check,indent=2)+'\n')
files={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()}
print(json.dumps(dict(hostname=socket.gethostname(),files=files)))
'''%(REMOTE,COMMITS)
    raw=subprocess.check_output(SSH+['/root/particle8_2_runs/env/bin/python -c '+shlex.quote(code)],text=True)
    inventory=json.loads(raw);(report/'REMOTE_FILE_INVENTORY.json').write_text(json.dumps(inventory,indent=2)+'\n')
    transport=shlex.join(SSH[:-1])
    subprocess.run(['rsync','-az','-e',transport,SSH[-1]+':'+REMOTE+'/results/ppt/',str(out)+'/'],check=True)
    mismatch=[p for p,h in inventory['files'].items() if not (out/p).is_file() or digest(out/p)!=h]
    result=dict(all_pass=not mismatch,checked_files=len(inventory['files']),mismatches=mismatch,
        remote_hostname=inventory['hostname'],remote_root=REMOTE+'/results/ppt',local_root=str(out),
        inventory_sha256=digest(report/'REMOTE_FILE_INVENTORY.json'))
    (report/'REMOTE_LOCAL_CHECKSUM_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    assert not mismatch,mismatch
    lock=json.loads((ROOT/'particle_3d/reports/particle8_2a/UPSTREAM_LOCK.json').read_text())
    changed=[p for p,h in lock['sha256'].items() if not Path(p).is_file() or digest(p)!=h]
    preserved=dict(all_pass=not changed,checked_files=len(lock['sha256']),changed=changed)
    (report/'UPSTREAM_IMMUTABILITY_VERIFICATION.json').write_text(json.dumps(preserved,indent=2)+'\n')
    assert not changed,changed
    print(json.dumps(dict(sync=result,upstream=preserved),indent=2))


if __name__=='__main__':main()

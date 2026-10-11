"""Pull only this new run and its four figures; never mirror/delete old data."""
from pathlib import Path
import json,subprocess,shlex,hashlib
C=json.loads(Path('/home/lzy/projects/.pvce_v2_current.json').read_text());root=Path(C['root']);base=Path(C['local']);SSH=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
script='''from pathlib import Path
import json
r=Path('/workspace/y_clot/project');bases=[r/'verification/lps_pvce_continuum_and_near_onset_audit_v2/20261011T142003Z',r/'results/lps_pvce_continuum_and_near_onset_audit_v2/20261011T142003Z'];items=[]
for b in bases:
 for p in b.rglob('*'):
  if p.is_file() and p.name!='body.data':items.append(str(p.relative_to(r)))
print(json.dumps(sorted(items)))'''
files=json.loads(subprocess.check_output(SSH+['python3 -'],input=script.encode()))
with subprocess.Popen(SSH+['cd /workspace/y_clot/project && '+shlex.join(['tar','-cf','-','--',*files])],stdout=subprocess.PIPE) as src:
 subprocess.run(['tar','-xf','-','-C',str(root)],stdin=src.stdout,check=True);src.stdout.close();assert src.wait()==0
receipt=[dict(path=f,bytes=(root/f).stat().st_size,sha256=hashlib.sha256((root/f).read_bytes()).hexdigest()) for f in files]
(base/'provenance/WSL_COMPACT_COLLECTION.json').write_text(json.dumps(dict(files=receipt,total_bytes=sum(r['bytes'] for r in receipt),excluded='Three text body.data files (full intact reference table) remain on server; compressed full reference arrays + exact writer + actual native inputs collected. No unrelated historical data copied.'),indent=2));print('collected',len(files),'files',sum(r['bytes'] for r in receipt),'bytes')

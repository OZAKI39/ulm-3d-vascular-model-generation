"""Copy only final documents/media/evidence to the AF server directory; verify bytes."""
import subprocess,shlex,time
from scripts.straight_tube_exp1.f10af.common import *
from scripts.straight_tube_exp1.f10af.remote import SSH,SERVER,call
def reseal():
    m=read(BASE/'INPUT_OUTPUT_SHA256.json')
    m['files']={str(p.relative_to(BASE)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*') if p.is_file() and p.name not in ['INPUT_OUTPUT_SHA256.json','LOCAL_RECEIPT.json'] and not p.name.endswith('.tar.gz') and (p.suffix!='.lock' or p.name.endswith('_environment.lock'))}
    m['identity']=identity();write(BASE/'INPUT_OUTPUT_SHA256.json',m)
start=time.perf_counter()
remote_base=SERVER+'/results/straight_tube_exp1/f10af/'
state=call('import subprocess,json\nr=subprocess.run(["supervisorctl","status"],stdout=subprocess.PIPE,text=True)\nprint(json.dumps([line for line in r.stdout.splitlines() if line.startswith("mevo_f10af_")]))')
import json
state=json.loads(state);assert all('RUNNING' not in s for s in state),state
write(BASE/'evidence/SUPERVISOR_FINAL_STATE.json',dict(rows=state,scope='Only F10AF programs; no unrelated job modified',identity=identity()))
selected=list(BASE.glob('*.md'))+[BASE/n for n in ['COMPUTE_COST.json','REPORT_COMPLETE.json','FINAL_GATE.json','INITIAL_STATE_PROVENANCE.json','REMOTE_RAW_TRAJECTORIES.json']]
selected += [p for name in ['figures','videos','evidence','logs'] for p in (BASE/name).rglob('*') if p.is_file() and not p.name.endswith('.tar.gz') and p.name not in ['REMOTE_DELIVERY_RECEIPT.json','REMOTE_DELIVERY_PLAN.json','REMOTE_DELIVERY_FILES.txt']]
files={str(p.relative_to(BASE)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(set(selected))}
write(BASE/'evidence/REMOTE_DELIVERY_PLAN.json',dict(files=files,server_root=remote_base,scope='Final documents, media and evidence only. No hf/, lf/, caches or historical results overwritten.',identity=identity()))
reseal()
names=list(files)+['evidence/REMOTE_DELIVERY_PLAN.json','INPUT_OUTPUT_SHA256.json'];listing=Path('/tmp/f10af_remote_delivery_files.txt');listing.write_text('\n'.join(names)+'\n')
subprocess.run(['rsync','-az','--files-from='+str(listing),'-e',shlex.join(SSH[:-1]),str(BASE)+'/',SSH[-1]+':'+remote_base],check=True)
code='root='+repr(remote_base)+'''\nfrom pathlib import Path
import json,hashlib
root=Path(root);plan=json.loads((root/'evidence/REMOTE_DELIVERY_PLAN.json').read_text());bad=[n for n,d in plan['files'].items() if hashlib.sha256((root/n).read_bytes()).hexdigest()!=d['sha256']]
assert not bad,bad
print(json.dumps(dict(PASS=True,checked=len(plan['files']),bytes=sum(d['bytes'] for d in plan['files'].values()))))
'''
receipt=json.loads(call(code));receipt.update(wall_s=time.perf_counter()-start,identity=identity());write(BASE/'evidence/REMOTE_DELIVERY_RECEIPT.json',receipt)
reseal()
listing.write_text('evidence/REMOTE_DELIVERY_RECEIPT.json\nINPUT_OUTPUT_SHA256.json\n')
subprocess.run(['rsync','-az','--files-from='+str(listing),'-e',shlex.join(SSH[:-1]),str(BASE)+'/',SSH[-1]+':'+remote_base],check=True)
print(receipt)

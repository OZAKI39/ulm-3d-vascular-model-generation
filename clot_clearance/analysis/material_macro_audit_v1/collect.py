"""Collect an explicit compact audit delivery; leave raw scalar histories in place."""
import hashlib,json,shlex,subprocess
from pathlib import Path

C=json.loads(Path('/home/lzy/projects/.macro_audit_current.json').read_text())
SSH=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
relative='verification/lammps_material_macro_calibration_and_timestep_audit_v1/'+C['run_id']
figrelative='results/lammps_material_macro_calibration_and_timestep_audit_v1/'+C['run_id']
script='''
import json,hashlib
from pathlib import Path
root=Path('/workspace/y_clot/project')
b=root/RELATIVE
f=root/FIGRELATIVE
assert (b/'FINAL_STATUS.json').is_file() and (f/'OPEN_RESULTS.html').is_file()
selected=[]
for p in sorted(b.rglob('*')):
 if not p.is_file():continue
 if p.name in ['history.csv','body.data','native_energy_density.dat','COMPACT_EXPORT.json','FINAL_METADATA_RECEIPT.json']:continue
 if p.suffix in ['.h5','.bin','.so','.pyc']:continue
 if p.suffix in ['.json','.npz','.csv','.md','.txt','.diff','.py','.cpp','.sh','.lammps','.log']:selected.append(p)
selected+=sorted(p for p in f.iterdir() if p.suffix in ['.png','.svg','.html'])
rows={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()) for p in selected}
print(json.dumps(rows))
'''.replace('FIGRELATIVE',repr(figrelative)).replace('RELATIVE',repr(relative))
rows=json.loads(subprocess.check_output(SSH+[shlex.join(['/workspace/y_clot/envs/clot/bin/python','-c',script])],text=True))
assert all(p.startswith(relative+'/') or p.startswith(figrelative+'/') for p in rows)
assert all('..' not in Path(p).parts for p in rows)
changed=[p for p,r in rows.items() if not (Path(C['root'])/p).is_file() or hashlib.file_digest((Path(C['root'])/p).open('rb'),'sha256').hexdigest()!=r['sha256']]
if changed:
 cmd=shlex.join(['tar','-cf','-','-C','/workspace/y_clot/project','--',*changed])
 with subprocess.Popen(SSH+[cmd],stdout=subprocess.PIPE) as remote:
  subprocess.run(['tar','-xf','-','-C',C['root']],stdin=remote.stdout,check=True)
  remote.stdout.close()
  assert remote.wait()==0
for p,r in rows.items():
 local=Path(C['root'])/p
 assert local.stat().st_size==r['bytes']
 assert hashlib.file_digest(local.open('rb'),'sha256').hexdigest()==r['sha256'],p
receipt=dict(status='COPIED_AND_SHA256_VERIFIED',server_root='/workspace/y_clot/project',local_root=C['root'],files=rows,total_bytes=sum(r['bytes'] for r in rows.values()),raw_data_retained_on_server=['*/history.csv','*/body.data','*/native_energy_density.dat','microbubble/dt_half/trajectory.h5'],note='review_curves.npz explicitly describes analysis sampling; all accepted scalar states remain in server history.csv. Later document/status metadata have a separate final receipt.')
(Path(C['local'])/'COMPACT_EXPORT.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(dict(files=len(rows),updated_files=len(changed),bytes=receipt['total_bytes'],status=receipt['status'])))

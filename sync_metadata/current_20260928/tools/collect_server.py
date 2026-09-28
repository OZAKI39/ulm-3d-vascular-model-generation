"""Read-only SSH inventory; collect current server evidence not already in snapshot."""
from pathlib import Path
import subprocess, json, csv, hashlib, collections, tarfile, os, sys
D=Path(__file__).resolve().parents[3]; M=D/'sync_metadata/current_20260928'
ROOTS={'flow_roi_only': '/workspace/flow_roi_only_balance_20260928', 'microbubble_roi_only': '/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928', 'flow_best_balance': '/workspace/flow_best_feasible_balance_20260928', 'microbubble_best_balance': '/workspace/microbubble_best_balance_dt0p5ms_n1500_20260928', 'rbc_work': '/workspace/hemocell_restore/work/rbc_stage1_20260915_165632', 'rbc_results': '/workspace/hemocell_restore/results/rbc_stage1_20260915_165632', 'solver_build': '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports/svmp_reuse_build.json', 'solver_environment': '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3n/use_gpu13_env.sh', 'solver_source': '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse'}

REMOTE=r'''
import os,sys,json,hashlib,stat
from pathlib import Path
roots=ROOTS_LITERAL
rows=[];omit=[];links=[];missing=[]
skip={'.git','.venv','env','__pycache__','.pytest_cache','build','bin','tmp','frozen_bundle','.cache'}
for label,root in roots.items():
 p=Path(root)
 if not p.exists():missing.append(root);continue
 def walk(p):
  if p.is_symlink():links.append(dict(source=str(p),target=os.readlink(p)));return
  if p.name in skip or p.name in ['USER_REQUEST.txt','supervisord.pid','.env','id_rsa','id_ed25519'] or p.suffix in ['.pyc','.o','.a','.so','.sock'] or p.name.endswith(('.tar.gz','.tar.zst','.zip')):
   omit.append(dict(source=str(p),reason='runtime, cache, compiled output, duplicate archive or historical request'));return
  if p.is_dir():
   for child in sorted(p.iterdir()):walk(child)
   return
  if not p.is_file():return
  n=p.stat().st_size
  if n>=50*1024**2 and not (p.name=='gpu_mesh_input.npz' and n<95*1024**2):omit.append(dict(source=str(p),bytes=n,reason='large file >=50 MiB'));return
  h=hashlib.sha256()
  with p.open('rb') as f:
   while b:=f.read(4*1024**2):h.update(b)
  rows.append(dict(label=label,source=str(p),relative=str(p.relative_to(Path(root))) if Path(root).is_dir() else p.name,bytes=n,sha256=h.hexdigest(),executable=bool(p.stat().st_mode&0o111)))
 walk(p)
binary=Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics')
binfo={}
if binary.exists():binfo=dict(path=str(binary),bytes=binary.stat().st_size,sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),included=False)
print(json.dumps(dict(files=rows,omitted=omit,links=links,missing_roots=missing,solver_binary=binfo)))
'''.replace('ROOTS_LITERAL',repr(ROOTS))

def main():
 result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','vast4090','python3 -'],input=REMOTE,text=True,capture_output=True,check=True,timeout=180)
 (M/'validation/server_inventory_ssh.log').write_text(result.stderr)
 inv=json.loads(result.stdout);(M/'server_inventory.json').write_text(json.dumps(inv,indent=2)+'\n')
 local=collections.defaultdict(list)
 for row in csv.DictReader((M/'local_file_inventory.csv').open()):local[row['sha256']].append(row['destination'])
 pointers={r['sha256']:r for r in json.loads((M/'upstream_lfs_pointers.json').read_text())}
 pointers.update({r['sha256']:r for r in json.loads((D/'sync_metadata/current_20260927/upstream_lfs_pointers.json').read_text())})
 mappings=[];selected={}
 for row in inv['files']:
  h=row['sha256'];entry=dict(row)
  if h in pointers:
   entry.update(disposition='upstream_lfs_pointer_metadata_only',snapshot=('sync_metadata/current_20260928/upstream_lfs_pointers.json' if any(r['sha256']==h for r in json.loads((M/'upstream_lfs_pointers.json').read_text())) else 'sync_metadata/current_20260927/upstream_lfs_pointers.json'),pointer_sha256=h)
  elif h in local:
   entry.update(disposition='identical_snapshot_blob',snapshot=local[h][0])
  elif h in selected:
   entry.update(disposition='identical_server_blob',snapshot=selected[h]['snapshot'])
  else:
   dest=Path('server_evidence/current_20260928')/row['label']/row['relative']
   # Server-local Git attributes must not affect subsequent add operations.
   if dest.name in ['.gitattributes','.gitignore','.ignore']:dest=dest.with_name(dest.name+'.archived')
   entry.update(disposition='collected',snapshot=str(dest));selected[h]=entry
  mappings.append(entry)
 (M/'server_file_map.json').write_text(json.dumps(mappings,indent=2)+'\n')
 summary=dict(remote_files=len(mappings),remote_bytes=sum(r['bytes'] for r in mappings),new_unique_files=len(selected),new_unique_bytes=sum(r['bytes'] for r in selected.values()),missing_roots=inv['missing_roots'])
 print(json.dumps(summary,indent=2),flush=True)
 archive=D.parent/'github_sync_20260928_server_evidence.tar'
 paths='\0'.join(r['source'].lstrip('/') for r in selected.values())+'\0'
 with archive.open('wb') as out:
  result=subprocess.run(['ssh','-o','BatchMode=yes','vast4090','tar -C / --null -T - -cf -'],input=paths.encode(),stdout=out,stderr=subprocess.PIPE,check=True)
 (M/'validation/server_collect_ssh.log').write_bytes(result.stderr)
 bysource={r['source'].lstrip('/'):r for r in selected.values()}
 with tarfile.open(archive) as tar:
  for member in tar:
   if not member.isfile():raise ValueError('unexpected tar member')
   row=bysource.pop(member.name);dst=D/row['snapshot'];dst.parent.mkdir(parents=True,exist_ok=True)
   payload=tar.extractfile(member).read()
   assert hashlib.sha256(payload).hexdigest()==row['sha256'],row['source']
   dst.write_bytes(payload);dst.chmod(0o755 if row['executable'] else 0o644)
 assert not bysource
 summary['collected_sha256_verified']=True
 (M/'server_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('Server collection verified.',flush=True)
 archive.unlink()

if __name__=='__main__':main()

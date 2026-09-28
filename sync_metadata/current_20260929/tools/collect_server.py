"""Read and copy only; never execute CFD, trajectories, or remote scripts."""
from pathlib import Path
import collections,hashlib,json,subprocess
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260929'
REMOTE='/workspace/brava_flow_roi_18mlmin_20260928'
code=r'''
from pathlib import Path
import hashlib,json,os
r=Path('/workspace/brava_flow_roi_18mlmin_20260928');rows=[];skipped=[]
for root,dirs,files in os.walk(r,followlinks=False):
 for name in dirs[:]:
  if name in {'.git','__pycache__','.cache','.venv'}:skipped.append(str(Path(root,name).relative_to(r)));dirs.remove(name)
 for name in files:
  p=Path(root,name);rel=str(p.relative_to(r))
  if p.is_symlink():rows.append(dict(path=rel,kind='symlink',target=os.readlink(p)));continue
  if not p.is_file() or p.suffix in {'.pid','.sock','.pyc'}:skipped.append(rel);continue
  with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
  size=p.stat().st_size
  pointer=p.read_text() if size<256 and p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1\n') else None
  reason=('Upstream LFS placeholder, payload absent' if pointer else 'Rebuildable compiled binary' if p.suffix in {'.so','.o','.a'} or name=='svmultiphysics' else 'Large reconstructible GPU input archive' if size>=95*1024**2 else '')
  row=dict(path=rel,kind='file',bytes=size,sha256=h,excluded_reason=reason)
  if pointer:row['lfs_pointer_text']=pointer
  rows.append(row)
print(json.dumps(dict(root=str(r),files=rows,skipped_runtime=skipped)))
'''
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 raw=subprocess.check_output(['ssh','-q','vast4090','/root/particle8_2_runs/env/bin/python -'],input=code,text=True)
 inv=json.loads(raw);(M/'server_inventory.json').write_text(json.dumps(inv,indent=2)+'\n')
 mapping=[];copy_rows=[]
 for row in inv['files']:
  if row['kind']=='symlink':
   mapping.append(dict(server_path=row['path'],action='LINK_RECORDED',target=row['target']));continue
  if row['excluded_reason']:
   mapping.append(dict(server_path=row['path'],action='EXCLUDED',reason=row['excluded_reason'],sha256=row['sha256']));continue
  rel=Path('brava_flow_roi_18mlmin')/row['path'];p=D/rel
  if p.name in {'.gitattributes','.gitignore','.ignore'}:
   rel=Path('sync_metadata/current_20260929/original_git_rules')/rel;p=D/rel
  if p.is_file() and sha(p)==row['sha256']:
   mapping.append(dict(server_path=row['path'],destination=str(rel),sha256=row['sha256'],action='MATCHED_LOCAL'));continue
  if p.exists():rel=Path('server_evidence/brava_current_20260929')/row['path']
  copy_rows.append((row,rel))
 # Transfer exact project-relative paths; rsync argument list avoids shell
 # interpolation of any logged content or authentication data.
 if copy_rows:
  temp=M/'validation/server_transfer'
  temp.mkdir(exist_ok=True)
  names='\n'.join(row['path'] for row,rel in copy_rows)+'\n'
  subprocess.run(['rsync','-az','-e','ssh -q','--files-from=-','vast4090:'+REMOTE+'/',str(temp)+'/'],input=names,text=True,check=True)
  import shutil
  for row,rel in copy_rows:
   src=temp/row['path'];assert sha(src)==row['sha256']
   dst=D/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.move(src,dst)
   mapping.append(dict(server_path=row['path'],destination=str(rel),sha256=row['sha256'],action='COPIED_SERVER'))
  shutil.rmtree(temp)
 (M/'server_file_map.json').write_text(json.dumps(mapping,indent=2)+'\n')
 result=subprocess.run(['ssh','-q','vast4090','supervisorctl -c '+REMOTE+'/supervisord.conf status brava_microbubble_sequence brava_flow_render'],text=True,capture_output=True)
 assert result.returncode in (0,3),result.stderr
 state=result.stdout
 assert 'brava_microbubble_sequence' in state and 'STOPPED' in state
 (M/'validation/server_status.txt').write_text(state)
 summary=dict(inventory_files=len(inv['files']),actions=dict(collections.Counter(r['action'] for r in mapping)),copied_bytes=sum(row['bytes'] for row,rel in copy_rows),scientific_jobs_started=0)
 (M/'server_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

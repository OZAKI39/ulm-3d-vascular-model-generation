from pathlib import Path
import subprocess,json,hashlib,tarfile,shutil
old=Path('/home/lzy/projects/.github-print-sync-20260925');work=old/'ulm_3D_vascular';audit=Path('/home/lzy/projects/vascular_printing/migration/20260926/retired_sync_workspace');audit.mkdir(parents=True)
def git(*args):return subprocess.check_output(['git','-C',str(work),*args])
head=git('rev-parse','HEAD').decode().strip();receipt=json.loads((old/'push_verification.json').read_text());assert head==receipt['commit']
status=git('status','--porcelain=v1','-z');(audit/'status_before.bin').write_bytes(status)
(audit/'uncommitted.patch').write_bytes(git('diff','--binary','HEAD'))
changed=[p for p in git('diff','--name-only','-z','HEAD').decode().split('\0') if p]
untracked=[p for p in git('ls-files','--others','--exclude-standard','-z').decode().split('\0') if p]
files=sorted(set(changed+untracked));inventory=[]
with tarfile.open(audit/'uncommitted_original_files.tar.gz','w:gz') as tar:
 for name in files:
  path=work/name
  if path.is_file():
   tar.add(path,arcname=name,recursive=False);inventory.append({'file':name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
with tarfile.open(audit/'uncommitted_original_files.tar.gz','r:gz') as tar:
 for r in inventory:assert hashlib.sha256(tar.extractfile(r['file']).read()).hexdigest()==r['sha256']
for p in old.iterdir():
 if p.is_file():shutil.copy2(p,audit/p.name)
# A byte inventory documents exactly which obsolete duplicate working files
# are removed; committed versions remain in the original repository history.
with (audit/'removed_worktree_files.jsonl').open('w') as stream:
 for p in work.rglob('*'):
  if p.is_file() and not p.is_symlink():
   stream.write(json.dumps({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})+'\n')
result=subprocess.run(['git','-C','/home/lzy/projects/ulm_3D_vascular','worktree','remove','--force',str(work)],capture_output=True,text=True)
(audit/'worktree_remove.log').write_text(result.stdout+result.stderr)
if result.returncode:raise RuntimeError(result.stderr)
assert not work.exists()
for p in old.iterdir():
 if p.is_file():
  assert hashlib.sha256(p.read_bytes()).hexdigest()==hashlib.sha256((audit/p.name).read_bytes()).hexdigest();p.unlink()
old.rmdir()
(audit/'retirement.json').write_text(json.dumps({'status':'REMOVED_OBSOLETE_SYNC_WORKTREE','commit':head,'commit_retained_in_original_git':True,'uncommitted_files_saved_and_verified':inventory,'git_branch_deleted':False,'remote_changed':False},indent=2))
print('Retired obsolete printing-sync worktree; uncommitted files preserved:',len(inventory))

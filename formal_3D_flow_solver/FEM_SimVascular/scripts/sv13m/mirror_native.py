import json,shlex,subprocess,tarfile,hashlib
from pathlib import Path
from remote import ROOT,REMOTE,ssh_prefix
ssh=ssh_prefix();R=ROOT/'reports/sv1_3m';O=ROOT/'outputs/sv1_3m/native';O.mkdir(exist_ok=True)
d=json.loads((R/'remote/native_artifacts.json').read_text());files=[]
for f in d['artifacts']:
 dst=O/Path(f['remote_path']).name
 subprocess.run(['rsync','-a','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+f['remote_path'],str(dst)],check=True)
 h=hashlib.sha256(dst.read_bytes()).hexdigest();assert h==f['sha256'];files.append(dict(f,local_path=str(dst.relative_to(ROOT))))
(R/'native_artifact_mirror.json').write_text(json.dumps(dict(status='PASS',files=files),indent=2)+'\n')
source=ROOT/'external/compat_cuda/petsc-3.19.6-cuda123-ghostfix';assert not source.exists()
with tarfile.open(O/'petsc-3.19.6.tar.gz') as t:
 temp=ROOT/'outputs/sv1_3m/local_source_extract';temp.mkdir();t.extractall(temp,filter='data');(temp/'petsc-3.19.6').rename(source)
subprocess.run(['patch','--batch','--fuzz=0','-p1','-i',str(ROOT/'patches/sv1_3m/petsc319_cuda_ghost_backport.patch')],cwd=source,check=True)
expected=json.loads((R/'remote/repair_03_source_integrity.json').read_text())['patched_files']
assert all(hashlib.sha256((source/n).read_bytes()).hexdigest()==h for n,h in expected.items())
(R/'WSL_source_mirror.json').write_text(json.dumps(dict(status='PASS',source=str(source.relative_to(ROOT)),files_verified=len(expected),native_compile_location='remote',WSL_is_source_of_truth=True),indent=2)+'\n')
print('Native artifacts mirrored and final PETSc source reconstructed in WSL with all file hashes verified.')

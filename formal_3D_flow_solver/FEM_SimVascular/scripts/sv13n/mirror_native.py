"""Mirror adopted executables, installs and exact build provenance into WSL."""
import hashlib,json,shlex,subprocess
from pathlib import Path
from remote import ROOT,REMOTE,ssh_prefix,fetch
ssh=ssh_prefix();R=ROOT/'reports/sv1_3n';O=ROOT/'outputs/sv1_3n/native';O.mkdir(exist_ok=True)
d=json.loads((R/'remote/native_artifacts.json').read_text());files=[]
for f in d['artifacts']:
 dst=O/Path(f['remote_path']).name
 subprocess.run(['rsync','-a','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+f['remote_path'],str(dst)],check=True)
 assert hashlib.sha256(dst.read_bytes()).hexdigest()==f['sha256'];files.append(dict(f,local_path=str(dst.relative_to(ROOT))))
(R/'native_artifact_mirror.json').write_text(json.dumps(dict(status='PASS',files=files),indent=2)+'\n')
for p in (R/'remote').glob('petsc_*_build.json'):
 d=json.loads(p.read_text())
 if d['status']=='PASS' and d['candidate_wrapper'].startswith(REMOTE+'/'):fetch(Path(d['candidate_wrapper']).name,ROOT/'scripts/sv13n'/Path(d['candidate_wrapper']).name)
S=ROOT/'external/petsc325/svMultiPhysics-compat';verified=[]
for p in (R/'remote').glob('svmp_*_build.json'):
 d=json.loads(p.read_text())
 if d.get('status')!='PASS':continue
 # Current adopted canonical builds only; initial failed/replaced attempts are retained separately.
 if p.stem not in ('svmp_cpu_build','svmp_gpu13_build','svmp_gpu123_build'):continue
 attempt=Path(d['build']).name;backend=p.stem.removeprefix('svmp_').removesuffix('_build')
 e=json.loads((R/'remote'/f'svmp_{backend}_{attempt}_source_expected.json').read_text())
 local_source=ROOT/'external/svMultiPhysics' if d.get('source_unmodified') else S
 mismatch=[n for n,h in e['files'].items() if hashlib.sha256((local_source/n).read_bytes()).hexdigest()!=h]
 assert not mismatch,('WSL adopted source differs from native binary',backend,mismatch)
 verified.append(dict(backend=backend,files=len(e['files'])))
(R/'WSL_source_mirror.json').write_text(json.dumps(dict(status='PASS',WSL_is_source_of_truth=True,adopted_solver_source=str(S.relative_to(ROOT)),verified=verified),indent=2)+'\n')

"""Hash-checked local mirror of stage-local native builds and source archives."""
import hashlib,json,shlex,subprocess
from remote import ROOT,ssh_prefix
R=ROOT/'reports/sv1_3p';d=json.loads((R/'remote/native_artifacts.json').read_text());ssh=ssh_prefix()
for f in d['files']:
 p=ROOT/f['local_path'];p.parent.mkdir(parents=True,exist_ok=True)
 if p.is_file() and not p.is_symlink() and hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256']:continue
 subprocess.run(['rsync','-aLz','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+f['remote_path'],str(p)],check=True)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256']
(R/'native_artifact_mirror.json').write_text(json.dumps(d,indent=2)+'\n');print('Native artifacts verified:',len(d['files']))

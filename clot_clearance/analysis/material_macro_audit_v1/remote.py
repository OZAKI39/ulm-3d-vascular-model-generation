"""Explicit task-scoped transport; no deletion, credentials, or whole-tree sync."""
import json,subprocess,sys,shlex
from pathlib import Path
C=json.loads(Path('/home/lzy/projects/.macro_audit_current.json').read_text())
SSH=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
if sys.argv[1]=='exec':
 raise SystemExit(subprocess.call(SSH+[sys.argv[2]]))
if sys.argv[1]=='python':
 raise SystemExit(subprocess.call(SSH+['/workspace/y_clot/envs/clot/bin/python -'],stdin=sys.stdin))
if sys.argv[1]=='sync':
 files=sys.argv[2:]
 for rel in files:
  src=Path(C['root'])/rel;dest=Path('/workspace/y_clot/project')/rel
  subprocess.run(SSH+['mkdir -p '+shlex.quote(str(dest.parent))],check=True)
  # New task paths only; existing files must be prior versions of this task.
  subprocess.run(['scp','-P','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes',str(src),'root@95.3.33.46:'+str(dest)],check=True)
 p=Path(C['local'])/'provenance/SYNC_EXPLICIT.json'
 old=json.loads(p.read_text()) if p.exists() else []
 import hashlib
 old.extend({'path':r,'sha256':hashlib.sha256((Path(C['root'])/r).read_bytes()).hexdigest()} for r in files)
 p.write_text(json.dumps(old,indent=2))

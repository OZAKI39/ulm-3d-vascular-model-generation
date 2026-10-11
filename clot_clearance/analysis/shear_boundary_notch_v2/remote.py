"""Task-scoped explicit synchronization; commands use the existing private environment."""
from pathlib import Path
import json,sys,subprocess,shlex,hashlib
C=json.loads(Path('/home/lzy/projects/.shear_notch_v2_current.json').read_text())
SSH=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
if sys.argv[1]=='python':raise SystemExit(subprocess.call(SSH+['source /workspace/y_clot/env.sh && cd /workspace/y_clot/project && python -'],stdin=sys.stdin))
if sys.argv[1]=='exec':raise SystemExit(subprocess.call(SSH+[sys.argv[2]]))
if sys.argv[1]=='sync':
 files=sys.argv[2:];r=Path(C['root']);assert all((r/p).is_file() for p in files)
 with subprocess.Popen(['tar','-cf','-','-C',str(r),'--',*files],stdout=subprocess.PIPE) as src:
  subprocess.run(SSH+[shlex.join(['tar','-xf','-','-C','/workspace/y_clot/project'])],stdin=src.stdout,check=True);src.stdout.close();assert src.wait()==0
 p=Path(C['local'])/'provenance/EXPLICIT_SYNC.json';a=json.loads(p.read_text()) if p.exists() else [];a.extend(dict(path=n,sha256=hashlib.sha256((r/n).read_bytes()).hexdigest()) for n in files);p.write_text(json.dumps(a,indent=2));print('explicitly synchronized',len(files),'files')

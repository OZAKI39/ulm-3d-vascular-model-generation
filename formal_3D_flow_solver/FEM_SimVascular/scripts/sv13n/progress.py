"""Bounded read-only runtime heartbeat; does not infer acceptance from partial logs."""
import json
from pathlib import Path
from remote import ROOT,REMOTE,python
code='''import json,re,time
from pathlib import Path
B=Path(%r);files=sorted((B/'logs').glob('*.log'),key=lambda p:p.stat().st_mtime)
rows=[]
for p in files[-4:]:
 text=p.read_text(errors='replace')
 solves=re.findall(r'^\\s*NS\\s+\\d+-.*$',text,re.M)
 if solves:rows.append(dict(log=p.name,last=solves[-1],linear_rows=len(solves),bytes=p.stat().st_size))
print(json.dumps(dict(timestamp=time.time(),recent_flow_logs=rows)))
'''%REMOTE
r=python(code);r.check_returncode();d=json.loads(r.stdout)
(ROOT/'reports/sv1_3n/progress.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,ensure_ascii=False))

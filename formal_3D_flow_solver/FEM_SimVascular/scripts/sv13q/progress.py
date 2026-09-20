"""Read-only compact native progress; never starts a solver."""
import json
from remote import python,REMOTE
code='''from pathlib import Path
import json,re,time
b=Path(BASE_PATH);result={}
for p in sorted((b/'logs').glob('R*.log'),key=lambda p:p.stat().st_mtime)[-2:]:
 if p.stem=='recovery_probe_run':continue
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-150000));s=f.read().decode(errors='replace')
 rows=re.findall(r'^\\s*NS\\s+\\d+-.*$',s,re.M)
 traces=re.findall(r'^SV13Q_(?:BEGIN|END).*$',s,re.M)
 result[p.stem]=dict(last_row=rows[-1] if rows else None,last_policy=traces[-1] if traces else None,age_s=time.time()-p.stat().st_mtime)
print(json.dumps(result,indent=2))
'''.replace('BASE_PATH',repr(REMOTE))
r=python(code);r.check_returncode();print(r.stdout)

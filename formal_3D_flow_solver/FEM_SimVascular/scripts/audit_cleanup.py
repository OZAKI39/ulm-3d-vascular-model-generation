#!/usr/bin/env python3
"""Check all project-maintained text, allowing only the two designated audit documents."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import write_json,now
path=ROOT/'reports/sv1/cleanup_summary.json'
audit=json.loads(path.read_text())
allowed={'reports/sv0/REPORT.md','reports/sv1/cleanup_summary.json'}
text_types={'.py','.yaml','.yml','.md','.json','.txt','.log','.xml','.toml','.cfg','.ini'}
paths=[]
for folder in ('src','scripts','tests','configs','inputs','reports','logs','outputs'):
    paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix.lower() in text_types)
paths.extend(p for p in ROOT.iterdir() if p.is_file() and (p.suffix in text_types or p.name.startswith('.git')))
paths.append(ROOT/'external/svMultiPhysics_commit.txt')
matches=[];scanned=[]
for p in sorted(set(paths)):
    relative=str(p.relative_to(ROOT))
    if relative in allowed: continue
    try: content=p.read_text().lower()
    except UnicodeError: continue
    hits=[term for term in audit['banned_active_terms'] if term in content]
    if hits: matches.append({'file':relative,'terms':hits})
    scanned.append(relative)
audit['active_text_scan']={'status':'FAIL' if matches else 'PASS','timestamp':now(),'scanned_files':len(scanned),'matches':matches,'paths':scanned}
audit['scan_scope']='All maintained text under source, scripts, tests, configs, inputs, reports, logs, outputs, plus root metadata and solver commit; two designated history/audit documents excepted. Immutable third-party distributions/source/builds and binary audit archives are not active owned source.'
audit['status']='FAIL' if matches else 'PASS'
write_json(path,audit)
print('Cleanup audit:',audit['status'],'files:',len(scanned),'matches:',matches)
raise SystemExit(bool(matches))

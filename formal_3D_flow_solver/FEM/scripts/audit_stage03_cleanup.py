#!/usr/bin/env python3
"""Read-only negative audit of all planned removals and active dependency paths."""
import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json

p=argparse.ArgumentParser();p.add_argument('--site',choices=('wsl','remote'),required=True)
p.add_argument('--output',default=None);a=p.parse_args()
plan=json.loads((ROOT/'inputs/stage03/stage018_cleanup_plan.json').read_text())
present=[f for f in plan['remove_files']+plan['remove_directories'] if (ROOT/f).exists() or (ROOT/f).is_symlink()]
pattern='|'.join(re.escape(n) for n in plan['forbidden_package_names'])
if shutil.which('rg'):
    scan=subprocess.run(['rg','-n','-i',pattern,'src','scripts','tests','configs','remote',
                         '-g','!*.pyc','-g','!source_manifest.json','-g','!connection.local.json'],cwd=ROOT,capture_output=True,text=True)
    assert scan.returncode in (0,1)
    matches=scan.stdout.splitlines();command=scan.args
else:
    matches=[];command=['Python equivalent text inventory; hidden environments and caches excluded']
    for category in ('src','scripts','tests','configs','remote'):
        for directory,children,files in os.walk(ROOT/category):
            children[:]=[d for d in children if not d.startswith('.') and d!='__pycache__']
            for name in files:
                if name in ('source_manifest.json','connection.local.json') or name.endswith('.pyc'):continue
                path=Path(directory)/name
                try:lines=path.read_text().splitlines()
                except UnicodeDecodeError:continue
                matches += [f'{path.relative_to(ROOT)}:{i+1}:{line}' for i,line in enumerate(lines) if re.search(pattern,line,re.I)]
assert not present and not matches,(present,matches)
assert all(sha256(ROOT/f)==v['sha256'] for f,v in plan['retained_shared_modules'].items())
report={'timestamp':timestamp(),'status':'PASS','site':a.site,'planned_paths_still_present':present,
        'active_source_scan_matches':matches,'active_source_scan_command':command,
        'stage018_test_modules_remaining':[str(p.relative_to(ROOT)) for p in (ROOT/'tests').glob('test_stage018_*.py')],
        'former_adapter_dependency_available':importlib.util.find_spec(plan['dependency_removal']['local_python_package']) is not None,
        'plan_sha256':sha256(ROOT/'inputs/stage03/stage018_cleanup_plan.json')}
assert not report['stage018_test_modules_remaining']
fingerprint={str(f.relative_to(ROOT)):sha256(f) for f in (ROOT/'remote/.env/conda-meta').glob('*.json')}
if (ROOT/'remote/.env/bin/python').is_file():fingerprint['remote/.env/bin/python']=sha256(ROOT/'remote/.env/bin/python')
digest=hashlib.sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()
original=json.loads((ROOT/f'outputs/stage03/cleanup/{a.site}.json').read_text())
report['fem_fingerprint_sha256']=digest
report['fem_environment_unchanged']=digest==original['fem_fingerprint_sha256']
assert report['fem_environment_unchanged']
write_json(ROOT/(a.output or f'outputs/stage03/cleanup/{a.site}_audit.json'),report)
print(json.dumps(report,indent=2))

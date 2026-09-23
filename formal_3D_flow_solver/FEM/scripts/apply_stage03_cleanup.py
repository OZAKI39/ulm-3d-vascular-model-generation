#!/usr/bin/env python3
"""Apply the user-authorized, enumerated deletion plan inside this project only."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
import shutil
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json

p=argparse.ArgumentParser()
p.add_argument('--site',choices=('wsl','remote'),required=True)
a=p.parse_args()
plan=json.loads((ROOT/'inputs/stage03/stage018_cleanup_plan.json').read_text())
before={str(f.relative_to(ROOT)):sha256(f) for f in (ROOT/'remote/.env/conda-meta').glob('*.json')}
if (ROOT/'remote/.env/bin/python').is_file():
    before['remote/.env/bin/python']=sha256(ROOT/'remote/.env/bin/python')
deleted=[]
for relative in plan['remove_files']+plan['remove_directories']:
    path=ROOT/relative
    assert Path(relative).parts[0] in ('src','scripts','tests','configs','remote','reports','outputs','logs','inputs','archive')
    assert '..' not in Path(relative).parts and not Path(relative).is_absolute()
    assert not any(parent.is_symlink() for parent in path.parents if parent!=ROOT and parent.is_relative_to(ROOT))
    if path.is_symlink() or path.is_file():
        path.unlink();deleted.append(relative)
    elif path.is_dir():
        shutil.rmtree(path);deleted.append(relative)
for relative in plan['remove_files']:
    path=ROOT/relative
    if path.suffix=='.py':
        for cached in (path.parent/'__pycache__').glob(path.stem+'.*.pyc'):
            cached.unlink()
assert all(not (ROOT/f).exists() for f in plan['remove_files']+plan['remove_directories'])
assert all(sha256(ROOT/f)==h for f,h in before.items())
for f,info in plan['retained_shared_modules'].items():
    assert sha256(ROOT/f)==info['sha256']
if a.site=='wsl':
    assert all(sha256(ROOT/f)==info['sha256'] for f,info in plan['retained_stage018'].items())
packages={d.metadata['Name'].lower():d.version for d in importlib.metadata.distributions()}
tool_names=set(plan['forbidden_package_names'])
conda_names=[]
for f in (ROOT/'remote/.env/conda-meta').glob('*.json'):
    conda_names.append(json.loads(f.read_text())['name'].lower())
assert not tool_names.intersection(packages) and not tool_names.intersection(conda_names)
out={'status':'PASS','site':a.site,'timestamp':timestamp(),'hostname':platform.node(),
     'plan_sha256':sha256(ROOT/'inputs/stage03/stage018_cleanup_plan.json'),'deleted':deleted,
     'all_planned_paths_absent':True,'production_environment_tool_dependencies_absent':True,
     'fem_environment_unchanged':True,'fem_fingerprint_count':len(before),
     'fem_fingerprint_sha256':hashlib.sha256(json.dumps(before,sort_keys=True).encode()).hexdigest(),
     'no_code_backups_created':True}
write_json(ROOT/f'outputs/stage03/cleanup/{a.site}.json',out)
print(json.dumps(out,indent=2))

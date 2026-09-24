#!/usr/bin/env python3
"""Repeat exact historical regression commands, then run the actual P6 suite."""
from pathlib import Path
import json,subprocess,sys
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle6_checkpoint import source_identity
from particle_3d.particle3_cases import write_json
LOGS=PACKAGE/'reports/particle6/logs'

def main():
    write_json(LOGS/'final_tested_sources.json',source_identity(REPO))
    commands=json.loads((LOGS/'baseline_commands.json').read_text());results=[]
    for row in commands:
        name=row['name'].replace('baseline_','final_')
        cmd=[x.replace('baseline_','final_') for x in row['command']]
        with (LOGS/(name+'.log')).open('w') as out:r=subprocess.run(cmd,cwd=row['cwd'],stdout=out,stderr=subprocess.STDOUT)
        results.append(dict(name=name,cwd=row['cwd'],command=cmd,returncode=r.returncode));write_json(LOGS/'final_commands.json',results)
        print(name,r.returncode,flush=True)
        if r.returncode:raise SystemExit(r.returncode)
    name='final_particle6';cmd=[str(REPO/'.venv/bin/python'),'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={LOGS/name}.xml','particle_3d/tests/particle6']
    with (LOGS/(name+'.log')).open('w') as out:r=subprocess.run(cmd,cwd=REPO,stdout=out,stderr=subprocess.STDOUT)
    results.append(dict(name=name,cwd=str(REPO),command=cmd,returncode=r.returncode));write_json(LOGS/'final_commands.json',results)
    print(name,r.returncode,flush=True)
    assert source_identity(REPO)['source_sha256']==json.loads((LOGS/'final_tested_sources.json').read_text())['source_sha256'],'Source changed during final tests'
    raise SystemExit(r.returncode)

if __name__=='__main__':main()

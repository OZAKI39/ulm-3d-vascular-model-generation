#!/usr/bin/env python3
"""Independent remote sync/run/fetch; cannot select the old FEM destination."""
import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import now,write_json
config=json.loads((ROOT/'remote/connection.local.json').read_text());destination=config['root']
assert destination.endswith('FEM_SimVascular_lzy') and destination!=config['old_fem_root']
ssh=['ssh','-p',str(config['port']),'-o','BatchMode=yes','-o','ConnectTimeout=15','-o','ServerAliveInterval=30','-o','ServerAliveCountMax=6','-o','IdentitiesOnly=yes','-i',config['identity_file']]
def run(args,logpath=None):
    p=subprocess.run(args,cwd=ROOT,text=True,capture_output=True)
    if logpath:
        write_json(logpath,{'timestamp':now(),'command':args,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
    print(p.stdout[-12000:]);print(p.stderr[-3000:],file=sys.stderr)
    if p.returncode:raise SystemExit(p.returncode)
parser=argparse.ArgumentParser();parser.add_argument('operation',choices=['init','sync','run','fetch']);parser.add_argument('--label',default='run');parser.add_argument('command',nargs=argparse.REMAINDER);a=parser.parse_args()
if a.operation=='init':
    code="from pathlib import Path; p=Path("+repr(destination)+"); p.mkdir(parents=True,exist_ok=True); print(p)"
    run(ssh+[config['host'],'python3 -B -c '+shlex.quote(code)],ROOT/'logs/sv0/remote_init.json')
elif a.operation=='sync':
    run(['rsync','-az','--relative','--exclude=__pycache__/','--exclude=.*','-e',shlex.join(ssh),*[str(ROOT)+'/./'+name for name in ['src','scripts','configs','inputs','tests','README.md','pyproject.toml']],config['host']+':'+destination+'/'])
elif a.operation=='run':
    args=a.command
    if args and args[0]=='--':args=args[1:]
    assert args
    remote_command='cd '+shlex.quote(destination)+' && PYTHONDONTWRITEBYTECODE=1 '+shlex.join(args)
    stamp=now().replace(':','').replace('.','_')
    run(ssh+[config['host'],remote_command],ROOT/f'logs/sv0/{stamp}_{a.label}.json')
elif a.operation=='fetch':
    target=ROOT/'outputs/sv0/remote_return';target.mkdir(parents=True,exist_ok=True)
    run(['rsync','-az','--include=/outputs/','--include=/outputs/sv0/***','--include=/reports/','--include=/reports/sv0/***','--include=/logs/','--include=/logs/sv0/***','--exclude=*','-e',shlex.join(ssh),config['host']+':'+destination+'/',str(target)+'/'])

#!/usr/bin/env python3
"""Read installed tools on each host, writing only inside this independent project."""
import argparse
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import now,write_json

def command(args,timeout=15):
    try:
        p=subprocess.run(args,text=True,capture_output=True,timeout=timeout)
        return {'args':args,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
    except (OSError,subprocess.TimeoutExpired) as exc:return {'args':args,'returncode':None,'error':str(exc)}

def probe():
    names=['simvascular','SimVascular','svpython','docker','podman','apptainer','singularity','skopeo','unshare','proot','git','git-lfs','python3']
    tools={n:shutil.which(n) for n in names}
    versions={n:command([p,'--version']) for n,p in tools.items() if p and n not in ('SimVascular','simvascular','svpython')}
    working={n:command([tools[n],'info']) for n in ('docker','podman') if tools[n]}
    candidates=[]
    for base in ('/usr/local','/opt','/usr/bin','/usr/local/bin'):
        p=Path(base)
        if p.exists():candidates += [str(x) for x in p.iterdir() if any(k in x.name.lower() for k in ('simvascular','svmultiphysics'))]
    spec=importlib.util.find_spec('sv')
    return {'timestamp':now(),'hostname':platform.node(),'architecture':platform.machine(),'os':platform.platform(),
            'os_release':Path('/etc/os-release').read_text() if Path('/etc/os-release').exists() else None,
            'python':{'executable':sys.executable,'version':sys.version},'executables':tools,'versions':versions,
            'container_runtime_info':working,'sv_module':None if spec is None else {'origin':spec.origin},
            'installation_candidates':candidates,'embedded_python':'NOT_FOUND' if not candidates else 'REQUIRES_API_CHECK',
            'cpu':command(['lscpu']),'memory':command(['free','-b']),'disk':command(['df','-B1',str(ROOT)]),
            'identity':command(['id']),'user_namespace_probe':command(['unshare','-Ur','true']) if tools['unshare'] else None,
            'fuse_device':Path('/dev/fuse').exists(),'docker_socket':Path('/var/run/docker.sock').exists(),
            'gpu_used':False,'old_fem_environment_used':False}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--site',choices=['wsl','remote'],required=True);a=p.parse_args()
    result=probe();result['site']=a.site
    write_json(ROOT/f'outputs/sv0/environment/{a.site}_probe.json',result)
    print(json.dumps({k:result[k] for k in ('site','hostname','architecture','executables','sv_module','installation_candidates','user_namespace_probe','container_runtime_info','fuse_device','docker_socket')},indent=2))

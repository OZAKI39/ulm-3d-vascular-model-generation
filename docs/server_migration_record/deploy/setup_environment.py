import json, os, subprocess, traceback
from pathlib import Path
from datetime import datetime, timezone

ROOT=Path('/workspace/brava_cfd')
steps=[]
def run(name,args,cwd=ROOT):
    print(name,flush=True)
    with (ROOT/'logs'/f'{name}.log').open('w') as f:
        p=subprocess.run(args,cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
    steps.append(dict(name=name,command=args,returncode=p.returncode))
    (ROOT/'reports/ENVIRONMENT_INSTALL_STEPS.json').write_text(json.dumps(steps,indent=2)+'\n')
    if p.returncode:raise RuntimeError(f'{name} failed: {p.returncode}')

result=dict(start_utc=datetime.now(timezone.utc).isoformat(),status='RUNNING')
try:
    os.environ['UV_CACHE_DIR']=str(ROOT/'runtime/uv-cache')
    run('create_cfd_venv',['uv','venv','--python','/usr/bin/python3','--seed',str(ROOT/'envs/cfd')])
    run('install_cfd_python',['uv','pip','install','--python',str(ROOT/'envs/cfd/bin/python'),'-r',str(ROOT/'runtime/requirements.in')])
    with (ROOT/'reports/requirements.lock.txt').open('w') as f:
        subprocess.run(['uv','pip','freeze','--python',str(ROOT/'envs/cfd/bin/python')],stdout=f,check=True)
    packages=ROOT/'deps/os-packages';packages.mkdir(exist_ok=True)
    run('download_isolated_mpi_dependency',['apt-get','download','libevent-pthreads-2.1-7t64'],packages)
    for package in packages.glob('*.deb'):
        run('extract_'+package.name,['dpkg-deb','-x',str(package),str(ROOT/'deps/sysroot')])
    result.update(status='PASS')
except Exception:
    result.update(status='FAIL',exception=traceback.format_exc())
    print(result['exception'],flush=True)
finally:
    result['end_utc']=datetime.now(timezone.utc).isoformat()
    (ROOT/'reports/ENVIRONMENT_INSTALL.json').write_text(json.dumps(result,indent=2)+'\n')
raise SystemExit(result['status']!='PASS')

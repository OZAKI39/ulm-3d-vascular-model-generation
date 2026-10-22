import hashlib,json,os,shutil,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path

root=Path('/workspace/brava_cfd')
def command(args):
    p=subprocess.run(args,capture_output=True,text=True,timeout=30)
    return dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
assert json.loads((root/'reports/ENVIRONMENT_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((root/'reports/PAYLOAD_INTEGRITY.json').read_text())['status']=='PASS'
cleanup=[]
for directory,dirs,files in os.walk(root,followlinks=False):
    if '.rsync-partial' in dirs:
        p=Path(directory)/'.rsync-partial'
        assert not p.is_symlink()
        cleanup.append(dict(path=str(p),reason='Task-created interrupted transfer scratch'))
        shutil.rmtree(p);dirs.remove('.rsync-partial')
cache=root/'runtime/uv-cache'
if cache.exists():
    assert not cache.is_symlink()
    cleanup.append(dict(path=str(cache),reason='Disposable package download/cache after successful environment install and pip dependency check'))
    shutil.rmtree(cache)
cfg=(root/'envs/cfd/pyvenv.cfg').read_text()
assert 'include-system-site-packages = false' in cfg
assert all('/venv/main' not in p for p in sys.path)
r=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),dedicated_root=str(root),python_prefix=sys.prefix,python_paths=sys.path,pyvenv_config=cfg,cleanup=cleanup)
r['gpu']=command(['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used,utilization.gpu,compute_cap','--format=csv,noheader'])
r['disk']=command(['df','-B1',str(root)])
r['project_size']=command(['du','-sb',str(root)])
r['pip_check']=command(['uv','pip','check','--python',str(root/'envs/cfd/bin/python')])
assert r['pip_check']['returncode']==0
r['supervisor']=command(['supervisorctl','-c',str(root/'runtime/supervisord.conf'),'status'])
r['resource_limits']={p:(Path('/sys/fs/cgroup')/p).read_text().strip() for p in ['cpu/cpu.cfs_quota_us','cpu/cpu.cfs_period_us','memory/memory.limit_in_bytes','memory/memory.usage_in_bytes','memory/memory.failcnt','memory/memory.oom_control']}
r['runtime_files']={}
for p in [root/'activate.sh',*sorted((root/'bin').iterdir()),root/'runtime/build_native.py',root/'runtime/setup_environment.py',root/'runtime/validate_environment.py',root/'runtime/verify_payload.py']:
    r['runtime_files'][str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
r['status']='PASS'
(root/'reports/FINAL_SERVER_STATE.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['status','timestamp_utc','gpu','project_size','resource_limits']},indent=2))

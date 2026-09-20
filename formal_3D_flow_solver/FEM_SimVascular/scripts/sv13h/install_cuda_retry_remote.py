"""Correct the new prefix parent permission, retaining the failed first attempt."""
import json,os,pwd,subprocess,time,shutil
from pathlib import Path
from environment_remote import BASE,env,digest,snapshot
R=BASE/'reports';L=BASE/'logs';old=json.loads((R/'cuda12_install.json').read_text())
assert old['status']=='FAIL' and not (R/'cuda12_install_attempt1.json').exists()
assert 'Unable to write to' in (L/'cuda12_installer_detail.log').read_text()
for rel in ['reports/cuda12_install.json','reports/post_install_environment.json','reports/cuda13_preservation.json','logs/cuda12_installer_detail.log','logs/cuda12_install.log']:
    p=BASE/rel;shutil.copyfile(p,p.with_name(p.stem+'_attempt1'+p.suffix))
source=json.loads((R/'cuda12_source_manifest.json').read_text());assert digest(source['archive'])==source['sha256']
prefix=Path(old['prefix']);assert not list(prefix.iterdir())
parent=prefix.parent;stat=parent.stat();user=pwd.getpwnam('nobody')
start=time.monotonic()
try:
    os.chown(parent,user.pw_uid,user.pw_gid)
    with (L/'cuda12_install_attempt2.log').open('x') as f:r=subprocess.run(old['command'],cwd=BASE/'tmp_cuda12_unprivileged',env=env(),stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=900)
finally:
    os.chown(parent,stat.st_uid,stat.st_gid)
    for base,dirs,files in os.walk(prefix,followlinks=False):
        os.chown(base,0,0)
        for name in files+dirs:os.chown(Path(base)/name,0,0,follow_symlinks=False)
detail=Path('/tmp/cuda-installer.log').read_text(errors='replace');assert str(prefix) in detail
(L/'cuda12_installer_detail.log').write_text(detail)
manifest=dict(old,exit_code=r.returncode,wall_time_s=time.monotonic()-start,status='PASS' if r.returncode==0 else 'FAIL',attempt=2,correction='New project external parent temporarily owned by installer uid; restored after installation',log='logs/cuda12_install_attempt2.log')
(R/'cuda12_install.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest),flush=True)
after=snapshot();(R/'post_install_environment.json').write_text(json.dumps(after,indent=2))
before=json.loads((R/'pre_install_environment.json').read_text())
checks={k:before[k]==after[k] for k in ['driver_files','cuda_default_link','cuda13_prefix','cuda13_files','ld_configuration']}
(R/'cuda13_preservation.json').write_text(json.dumps({'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks},indent=2))
assert all(checks.values()),'SYSTEM_CUDA_OR_DRIVER_CHANGED'
assert r.returncode==0 and (prefix/'bin/nvcc').is_file(),'CUDA12_INSTALL_FAILED'
print('Toolkit-only custom-prefix installation PASS; existing CUDA13 and driver unchanged',flush=True)

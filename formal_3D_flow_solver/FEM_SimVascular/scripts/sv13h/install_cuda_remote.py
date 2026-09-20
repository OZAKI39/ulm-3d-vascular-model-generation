"""Run NVIDIA's toolkit-only installer without permission to alter system CUDA."""
import json,os,pwd,subprocess,time
from pathlib import Path
from environment_remote import BASE,env,digest,snapshot
R=BASE/'reports';L=BASE/'logs'
source=json.loads((R/'cuda12_source_manifest.json').read_text())
assert digest(source['archive'])==source['sha256']
prefix=BASE/'external/cuda-12.6';tmp=BASE/'tmp_cuda12_unprivileged'
assert not prefix.exists();assert not (R/'cuda12_install.json').exists()
user=pwd.getpwnam('nobody')
for p in (prefix,tmp):p.mkdir();os.chown(p,user.pw_uid,user.pw_gid)
# Only these new directories are writable to the installer identity.
# No user account, HOME, shell startup, system linker or CUDA symlink is edited.
cmd=['/usr/bin/setpriv','--reuid='+str(user.pw_uid),'--regid='+str(user.pw_gid),'--clear-groups','--no-new-privs','sh',source['archive'],'--silent','--toolkit','--toolkitpath='+str(prefix),'--defaultroot='+str(prefix),'--no-man-page','--tmpdir='+str(tmp)]
start=time.monotonic()
with (L/'cuda12_install.log').open('x') as f:r=subprocess.run(cmd,cwd=tmp,env=env(),stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=900)
manifest={'command':cmd,'exit_code':r.returncode,'wall_time_s':time.monotonic()-start,'prefix':str(prefix),'toolkit_only':True,'driver_install':False,'installer_uid':user.pw_uid,'source_sha256':source['sha256'],'status':'PASS' if r.returncode==0 else 'FAIL'}
(R/'cuda12_install.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest),flush=True)
for p in (tmp/'cuda-installer.log',Path('/tmp/cuda-installer.log')):
    if p.exists():
        text=p.read_text(errors='replace')
        if str(prefix) in text:(L/'cuda12_installer_detail.log').write_text(text)
after=snapshot();(R/'post_install_environment.json').write_text(json.dumps(after,indent=2))
before=json.loads((R/'pre_install_environment.json').read_text())
keys=['driver_files','cuda_default_link','cuda13_prefix','cuda13_files','ld_configuration']
checks={k:before[k]==after[k] for k in keys}
(R/'cuda13_preservation.json').write_text(json.dumps({'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks},indent=2))
assert all(checks.values()),'SYSTEM_CUDA_OR_DRIVER_CHANGED'
assert r.returncode==0,'CUDA12_INSTALL_FAILED'
assert (prefix/'bin/nvcc').is_file()
print('Toolkit installed; complete CUDA13/driver fingerprints unchanged',flush=True)

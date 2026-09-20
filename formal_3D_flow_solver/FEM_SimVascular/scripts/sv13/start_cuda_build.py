#!/usr/bin/env python3
import json,shlex,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from remote import ssh_prefix,REMOTE_ROOT,run_command
from sv_validation.sv13 import REPORT,LOG,load
from sv_validation.provenance import write_json,sha256,now
assert load('cpu_validation')['status']=='PASS' and load('gpu_environment')['status']=='AVAILABLE'
prefix=ssh_prefix();ssh=shlex.join(prefix[:-1]);host=prefix[-1]
files=[ROOT/'external/downloads/sv1_1/petsc-3.19.6.tar.gz',ROOT/'scripts/sv13/build_cuda_remote.py']
for file in files:
    p=subprocess.run(['rsync','-a','-e',ssh,str(file),host+':'+REMOTE_ROOT+'/'],capture_output=True,text=True)
    assert p.returncode==0,p.stderr
write_json(REPORT/'cuda_transfer.json',{'timestamp':now(),'files':[{'source':str(f.relative_to(ROOT)),'sha256':sha256(f),'destination':REMOTE_ROOT+'/'+f.name} for f in files]})
with (LOG/'cuda_remote_driver.log').open('x') as out:
    p=run_command(['/usr/bin/python3','-u','-B',REMOTE_ROOT+'/build_cuda_remote.py'],stdout=out,stderr=subprocess.STDOUT)
print('Remote CUDA build exit:',p.returncode,flush=True)
# Fetch all build logs even on failure; the local workspace owns the evidence.
dest=ROOT/'outputs/sv1_3/remote_build';dest.mkdir(parents=True,exist_ok=True)
for name in ('logs/','cuda_build_manifest.json','petsc-3.19.6/configure.log'):
    target=dest/('logs' if name=='logs/' else Path(name).name)
    if name=='logs/':target.mkdir(exist_ok=True)
    fetch=subprocess.run(['rsync','-a','-e',ssh,host+':'+REMOTE_ROOT+'/'+name,str(target)+('/' if name=='logs/' else '')],capture_output=True,text=True)
    if fetch.returncode:print('Fetch issue:',name,fetch.stderr,flush=True)
manifest=json.loads((dest/'cuda_build_manifest.json').read_text())
manifest['local_evidence_directory']=str(dest.relative_to(ROOT));manifest['remote_driver_exit_code']=p.returncode
write_json(REPORT/'cuda_petsc_build.json',manifest)
print(json.dumps({k:v for k,v in manifest.items() if k not in ('steps','configure_options','linked_libraries')},indent=2),flush=True)

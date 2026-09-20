"""Download the official runfile, verify NVIDIA's published MD5, freeze SHA256."""
import hashlib,json,subprocess,time,urllib.request
from pathlib import Path
from environment_remote import BASE,env
R=BASE/'reports';L=BASE/'logs';E=BASE/'external'
for p in (L,E):p.mkdir(exist_ok=True)
assert (R/'pre_install_environment.json').exists()
selection=json.loads((BASE/'configs/cuda_selection.json').read_text())
url=selection['url'];assert url.startswith('https://developer.download.nvidia.com/compute/cuda/12.6.3/local_installers/')
archive=E/url.rsplit('/',1)[1];assert not archive.exists()
md5=hashlib.md5();sha=hashlib.sha256();start=time.monotonic();size=0;next_notice=512*1024*1024
with urllib.request.urlopen(url,timeout=60) as src,archive.with_suffix('.part').open('wb') as f:
    length=int(src.headers['Content-Length'])
    for block in iter(lambda:src.read(8*1024*1024),b''):
        f.write(block);md5.update(block);sha.update(block);size+=len(block)
        if size>=next_notice:print(f'Downloaded {size/1024**3:.2f}/{length/1024**3:.2f} GiB',flush=True);next_notice+=512*1024*1024
assert size==length and md5.hexdigest()==selection['official_md5']
archive.with_suffix('.part').rename(archive)
manifest=dict(selection,sha256=sha.hexdigest(),actual_md5=md5.hexdigest(),bytes=size,wall_time_s=time.monotonic()-start,archive=str(archive),status='PASS',frozen_before_installer_execution=True,integrity='Official published MD5 matched; SHA256 computed and frozen before use (not claimed to be an NVIDIA-published SHA256)')
(R/'cuda12_source_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest),flush=True)
with (L/'cuda12_installer_help.log').open('w') as f:
    r=subprocess.run(['sh',str(archive),'--help'],env=env(),stdout=f,stderr=subprocess.STDOUT,timeout=60)
assert r.returncode==0
print('Installer help recorded; installation has not started',flush=True)

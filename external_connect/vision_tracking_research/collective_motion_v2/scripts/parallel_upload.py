"""Bounded transfer of independent archive chunks; no remote overwrite of sources."""
import argparse,concurrent.futures,subprocess,time,json,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--workers',type=int,default=8);a=p.parse_args();folder=a.archive.with_suffix('.parts');folder.mkdir(exist_ok=True)
with a.archive.open('rb') as f:
 i=0
 while b:=f.read(8*2**20):
  (folder/f'part_{i:04}').write_bytes(b);i+=1
ssh=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
remote='/workspace/f10v_vision_20261017/setup/'+folder.name
subprocess.run(ssh+['mkdir -p '+remote],check=True)
def upload(f):
 t=time.perf_counter();subprocess.run(['scp','-q','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283',str(f),'root@95.3.33.46:'+remote+'/'],check=True);print(f.name,time.perf_counter()-t,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(upload,sorted(folder.glob('part_*'))))
subprocess.run(ssh+[f'cat {remote}/part_* > /workspace/f10v_vision_20261017/setup/{a.archive.name}.complete'],check=True)
print('complete',hashlib.sha256(a.archive.read_bytes()).hexdigest(),flush=True)

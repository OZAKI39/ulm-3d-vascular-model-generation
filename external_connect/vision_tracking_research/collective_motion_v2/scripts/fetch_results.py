"""Fetch this stage's small artifacts and best/last weights with parallel chunks."""
import subprocess,concurrent.futures,hashlib,tarfile,tempfile,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ssh=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
base='/workspace/f10v_vision_20261017';archive=base+'/setup/delivery_core.tar.gz'
cmd=f"test -f {archive} || tar --exclude='epoch_*.pt' -czf {archive} -C {base}/runs/prototype_7aed evidence results training; mkdir -p {base}/setup/delivery_parts; split -b 8M -d -a 4 {archive} {base}/setup/delivery_parts/part_; sha256sum {archive}; ls {base}/setup/delivery_parts/part_*"
out=subprocess.check_output(ssh+[cmd],text=True);lines=out.splitlines();hashval=next(s.split()[0] for s in lines if 'delivery_core.tar.gz' in s);parts=[s for s in lines if s.startswith(base+'/setup/delivery_parts/part_')];dest=Path('/tmp/f10v_delivery_core');dest.mkdir(exist_ok=True)
def get(p):
 q=dest/Path(p).name
 for attempt in range(5):
  result=subprocess.run(['scp','-q','-o','ConnectTimeout=15','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283','root@95.3.33.46:'+p,str(q)])
  if result.returncode==0:return q
  time.sleep(2+attempt)
 raise RuntimeError('transfer failed '+p)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:files=list(pool.map(get,parts))
file=dest/'delivery_core.tar.gz'
with file.open('wb') as f:
 for p in sorted(files):f.write(p.read_bytes())
assert hashlib.sha256(file.read_bytes()).hexdigest()==hashval
with tarfile.open(file) as t:t.extractall(ROOT,filter='data')
(ROOT/'evidence/result_transfer.json').write_text(json.dumps(dict(archive_sha256=hashval,remote_archive=archive,parts=len(parts),excluded='periodic checkpoints retained on server; best and last copied'),indent=2))
print('verified/extracted',file,flush=True)

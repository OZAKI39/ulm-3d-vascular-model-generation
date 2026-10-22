from pathlib import Path
import subprocess,concurrent.futures,time
root=Path(__file__).resolve().parents[1];folder=root/'results/synthetic_demo';folder.mkdir(exist_ok=True);base='root@95.3.33.46:/workspace/f10v_vision_20261017/runs/prototype_7aed/results/synthetic_demo/'
files=[c+e for c in ['translation','rotation_cw','contraction','fixed_background','defocus','dropout','step'] for e in ['.npz','.json']]+['manifest.json']
def fetch(name):
 for attempt in range(5):
  r=subprocess.run(['scp','-q','-o','ConnectTimeout=15','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283',base+name,str(folder/name)])
  if r.returncode==0:return
  time.sleep(2)
 raise RuntimeError(name)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as p:list(p.map(fetch,files))

import sys,time,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import torch,numpy as np
from networks import build_model
from learning import CachedTriples,loss_fn
from vcommon import OUT,dump,identity,stats
from torch.utils.data import DataLoader

torch.set_num_threads(4);torch.backends.cuda.matmul.allow_tf32=True
rows=[];ds=CachedTriples('train')
for name in ['direct_fused','geo_fused']:
 for batch_size in [32,96,192]:
  b=next(iter(DataLoader(ds,batch_size=batch_size)))
  b={k:v.cuda() if torch.is_tensor(v) else v for k,v in b.items()}
  for amp in [False,True]:
   torch.manual_seed(42);m=build_model(name).cuda();opt=torch.optim.AdamW(m.parameters(),lr=.001);times=[];torch.cuda.reset_peak_memory_stats()
   for step in range(8):
    torch.cuda.synchronize();t=time.perf_counter();opt.zero_grad(set_to_none=True)
    with torch.autocast('cuda',dtype=torch.bfloat16,enabled=amp):out=m(b['images'],b['dt']);loss,parts=loss_fn(out,b)
    loss.backward();grad=torch.nn.utils.clip_grad_norm_(m.parameters(),10);opt.step();torch.cuda.synchronize()
    if step>=3:times.append(time.perf_counter()-t)
    assert torch.isfinite(loss) and torch.isfinite(grad)
   row=dict(name=name,batch=batch_size,amp='BF16' if amp else 'FP32',step_s=stats(times),samples_per_s=batch_size/np.mean(times),peak_mb=torch.cuda.max_memory_allocated()/2**20,last_loss=float(loss),gradient_norm=float(grad));rows.append(row);print(json.dumps(row),flush=True)
   del m,opt,out;torch.cuda.empty_cache()
# Unit-level parity uses untrained same weights, not accuracy/selection evidence.
dump(OUT/'evidence/smoke_profile.json',dict(identity=identity(),rows=rows,scope='development hardware/gradient check, not formal training'))

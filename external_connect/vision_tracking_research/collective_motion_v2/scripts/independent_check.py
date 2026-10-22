import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,torch
from torch.utils.data import DataLoader
from inference import Predictor
from learning import CachedTriples,loss_fn
from networks import build_model
from vcommon import ROOT,OUT,LIMITS,dump,sha,identity

def main():
    torch.set_num_threads(4);frozen=json.loads((OUT/'results/model_freeze.json').read_text());b=next(iter(DataLoader(CachedTriples('validation'),batch_size=128)));rows=[]
    for name in ['direct_gray','direct_fused','geo_fused']:
        key=frozen['champions_by_validation_only'][name];path=OUT/frozen['models'][key]['checkpoint'];m=Predictor(path);out=m.outputs(m.batch(b['images'],b['dt']),b['dt'].numpy());reference=np.load(OUT/'results'/f'{key}_validation.npz');diff=np.abs(out['mean']-reference['mean'][:128]);rows.append(dict(model=key,checkpoint_sha256=sha(path),max_absolute_difference=diff.max(0).tolist(),max_difference_over_tolerance=float(np.max(diff/LIMITS)),all_finite=bool(np.isfinite(out['mean']).all()),scope='independent local Torch 2.11 / RTX4060 versus server Torch2.14 / RTX4090; fixed validation batch; no retuning'))
        assert np.max(diff/LIMITS)<.1
    # Restore a last checkpoint in memory and verify one optimizer step; do not save or select it.
    path=OUT/'training/direct_fused_seed1729/last.pt';before=sha(path);c=torch.load(path,map_location='cpu',weights_only=False);m=build_model('direct_fused').cuda();m.load_state_dict(c['model']);opt=torch.optim.AdamW(m.parameters(),lr=.001);opt.load_state_dict(c['optimizer'])
    scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(opt,40);scheduler.load_state_dict(c['scheduler']);b=next(iter(DataLoader(CachedTriples('train'),batch_size=16)));b={k:v.cuda() for k,v in b.items()};p0=next(m.parameters()).detach().clone();m.train();opt.zero_grad();out=m(b['images'],b['dt']);loss,_=loss_fn(out,b);loss.backward();opt.step();changed=not torch.equal(p0,next(m.parameters()).detach());assert changed and sha(path)==before
    dump(OUT/'evidence/independent_checkpoint_check.json',dict(identity=identity(),models=rows,optimizer_resume_in_memory=dict(pass_check=changed,restored_epoch=c['epoch']+1,loss=float(loss.detach()),formal_file_unchanged=True,training_result_used=False),cuda=torch.version.cuda,torch=torch.__version__,device=torch.cuda.get_device_name(0)))
if __name__=='__main__':main()

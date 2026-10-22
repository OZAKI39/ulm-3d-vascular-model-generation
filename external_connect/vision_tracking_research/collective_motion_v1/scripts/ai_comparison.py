"""One official SEA-RAFT-S checkpoint, 24 endpoint pairs, no optimization/training."""
import sys,json,time,argparse,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
import cv2
from affine import fit_affine,parameters
from methods import preprocess

def run():
    import torch
    from types import SimpleNamespace
    repo=ROOT/'third_party/SEA-RAFT';sys.path.insert(0,str(repo/'core'))
    from raft import RAFT
    cfg=json.loads((repo/'config/eval/spring-S.json').read_text());args=SimpleNamespace(**cfg)
    torch.set_num_threads(1);torch.manual_seed(43)
    model_id='MemorySlices/Tartan-C-T-TSKH-spring540x960-S'
    from huggingface_hub import snapshot_download
    snapshot=Path(snapshot_download(model_id,local_dir=ROOT/'models/SEA-RAFT-S',allow_patterns=['*.json','*.safetensors','README.md']))
    model=RAFT.from_pretrained(str(snapshot),args=args).eval().cuda()
    def flow(a,b):
        aa=torch.from_numpy(np.repeat(a[None,None],3,axis=1).copy()).float().cuda()
        bb=torch.from_numpy(np.repeat(b[None,None],3,axis=1).copy()).float().cuda()
        with torch.inference_mode():ff=model(aa,bb,iters=4,test_mode=True)['flow'][-1]
        return ff[0].permute(1,2,0).cpu().numpy()
    cases=['translation','rotation_cw','contraction','shear','fixed_background','dense'];pairs=[]
    for name in cases:
        z=np.load(ROOT/f'results/synthetic/test/7201_{name}.npz');ims=z['frames'];HH=z['H']
        for window in P['windows']:
            end=15;begin=end-window+1;dt=(window-1)*P['synthetic_dt_s'];H=HH[end]@np.linalg.inv(HH[begin]);truth=parameters(H[:2,:2],H[:2,2],dt,np.array([128.,128.]))
            pairs.append({'source':'SYNTHETIC_MOTION_TEST','case':name,'seed':7201,'window':window,'begin':begin,'end':end,'truth':truth.tolist(),'dt_s':dt,'a':ims[begin],'b':ims[end]})
    z=np.load(ROOT/'results/real_frames.npz');ims=z['images'];pts=z['pts']
    for end in [20,90]:
        for window in P['windows']:
            begin=end-window+1;pairs.append({'source':'RECORDED_CHRONOS','case':'real','window':window,'begin':begin,'end':end,'truth':None,'dt_s':float(pts[end]-pts[begin]),'a':cv2.cvtColor(ims[begin],cv2.COLOR_BGR2GRAY),'b':cv2.cvtColor(ims[end],cv2.COLOR_BGR2GRAY)})
    # Limited warmup on first pair, explicitly excluded from reported steady-state samples.
    t=time.perf_counter();flow(pairs[0]['a'],pairs[0]['b']);torch.cuda.synchronize();warmup=time.perf_counter()-t
    torch.cuda.reset_peak_memory_stats();rows=[];start=time.perf_counter();resources=ResourceMonitor()
    for pair in pairs:
        a=pair.pop('a');b=pair.pop('b');h,w=a.shape
        torch.cuda.synchronize();tick=time.perf_counter();ff=flow(a,b);bb=flow(b,a);torch.cuda.synchronize()
        x,y=np.meshgrid(np.arange(16,w-15,16),np.arange(16,h-15,16));xy=np.c_[x.ravel(),y.ravel()].astype(float);d=ff[y,x]
        bx=(x+d[:,:,0]).astype(np.float32);by=(y+d[:,:,1]).astype(np.float32)
        rev=cv2.remap(bb,bx,by,cv2.INTER_LINEAR,borderMode=cv2.BORDER_CONSTANT,borderValue=float('nan'))
        fb=np.linalg.norm(d+rev,axis=2);signal=preprocess(a,'local');std=np.array([np.std(signal[yy-16:yy+16,xx-16:xx+16]) for xx,yy in zip(x.ravel(),y.ravel())]).reshape(x.shape)
        good=(fb<1.5)&(std>=P['piv']['texture_std_min'])&np.isfinite(d).all(2)
        algo=(time.perf_counter()-tick)*1000
        estimate=fit_affine(xy[good.ravel()],d[good],pair['dt_s'],a.shape,'SEA_RAFT_S',bootstrap=P['bootstrap_replicates'])
        row={**pair,'method':'SEA_RAFT_S','estimate':estimate,'algorithm_ms':algo,'valid_vector_count':int(good.sum()),'forward_backward_residual_median_px':float(np.nanmedian(fb)),'rejected_fraction':float(1-good.mean())};rows.append(row)
        np.savez_compressed(ROOT/f'results/ai_flow_{len(rows):02d}.npz',flow=ff,reverse=bb,valid=good)
        dump(ROOT/'results/ai_comparison.json',{'status':'RUNNING','rows':rows})
        print('AI',len(rows),pair['case'],pair['window'],round(algo,2),estimate['observability_status'],flush=True)
        if time.perf_counter()-start>300:raise RuntimeError('AI bounded 5 minute inference budget')
    resource_data=resources.close()
    dump(ROOT/'results/ai_comparison.json',{'status':'COMPLETE','identity':identity(),'resources':resource_data,'official_git':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),'model_id':model_id,
      'model_files_sha256':{str(p.relative_to(ROOT)):sha(p) for p in snapshot.glob('*') if p.is_file()},'code_license':'BSD-3-Clause','weight_license':'BSD-3-Clause per official model card','torch':torch.__version__,'cuda':torch.version.cuda,'device':torch.cuda.get_device_name(),'input':'native 256x256 / 256x384 grayscale repeated three channels, 0..255; official internal normalization and padding; no input resize','iterations':4,
      'warmup_s':warmup,'elapsed_s':time.perf_counter()-start,'allocated_peak_mb':torch.cuda.max_memory_allocated()/2**20,'reserved_peak_mb':torch.cuda.max_memory_reserved()/2**20,'rows':rows,'timing_ms':stats([r['algorithm_ms'] for r in rows]),'pair_count':len(rows),'model_forward_calls':1+2*len(rows),'fine_tuning':False,'admission':'NOT_ADMITTED_PENDING_TRUTH_COMPARISON'})

if __name__=='__main__':run()

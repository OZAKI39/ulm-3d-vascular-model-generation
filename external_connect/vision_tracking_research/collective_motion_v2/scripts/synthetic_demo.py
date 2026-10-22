"""Supplemental continuous demo from unmodified F10U generator, frozen models."""
import sys,os,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,torch
from legacy import read_frames,original_synthetic,piv
from heatmaps import HeatmapDetector,centers_to_heatmap
from inference import Predictor,estimate_from_prediction
from geometry import endpoint_six
from vcommon import ROOT,OUT,P,dump,sha,identity

def main():
    torch.set_num_threads(4);freeze=json.loads((OUT/'results/model_freeze.json').read_text());protocol=json.loads((ROOT/'configs/demo_protocol.json').read_text());video=os.environ.get('F10V_VIDEO',P['source_video'])
    source=read_frames([protocol['source_frame']],video)[protocol['source_frame']];models={k:Predictor(OUT/freeze['models'][freeze['champions_by_validation_only'][k]]['checkpoint']) for k in ['direct_gray','direct_fused','geo_fused']}
    detector=HeatmapDetector();folder=OUT/'results/synthetic_demo';folder.mkdir(parents=True,exist_ok=True);start=time.perf_counter();manifest=[]
    for name,seed in zip(protocol['cases'],protocol['seeds']):
        frames,meta=original_synthetic.build_sequence(name,seed,source,length=protocol['frames_per_case']);records=[];objects=[];heats=[];badheats=[]
        for k,im in enumerate(frames):
            o,g,h=detector(im);objects.append(o);heats.append(h);rng=np.random.default_rng(seed+k+700000);xy=np.array([z['xy'] for z in o]).reshape(-1,2);xy=xy[rng.random(len(xy))>.6];xy=np.r_[xy,rng.uniform(0,256,(20,2))];badheats.append(centers_to_heatmap(xy,(256,256)))
            row=dict(frame=k,time_s=k*meta['dt'],objects=o,truth=None,predictions={})
            if k>=2:
                dt=2*meta['dt'];H=meta['H'][k]@np.linalg.inv(meta['H'][k-2]);row['truth']=endpoint_six(H,dt,im.shape).tolist();ts=np.arange(k-2,k+1)*meta['dt']
                r=piv(frames[k-2],im,dt,mode='candidate',objects_a=objects[k-2],objects_b=o);row['predictions']['PIV_SINGLE']=r['estimate'];row['piv_vectors']=r['vectors']
                for key,m in models.items():
                    pred,timing=m.frames(frames[k-2:k+1],heats[k-2:k+1],ts);cal=freeze['models'][freeze['champions_by_validation_only'][key]];row['predictions'][key]=estimate_from_prediction(pred,cal,key,m.sha)
                    if key=='direct_fused':
                        pred,_=m.frames(frames[k-2:k+1],badheats[k-2:k+1],ts);row['predictions']['direct_fused_corrupt']=estimate_from_prediction(pred,cal,'direct_fused_corrupt',m.sha)
            records.append(row)
        np.savez_compressed(folder/f'{name}.npz',frames=frames,H=meta['H'],heat=np.uint8(np.clip(np.array(heats)*255,0,255)),corrupt_heat=np.uint8(np.clip(np.array(badheats)*255,0,255)))
        dump(folder/f'{name}.json',dict(case=name,seed=seed,source_frame=protocol['source_frame'],source_pts=source[1],dt_s=meta['dt'],scope=protocol['generator'],records=records,identity=identity()))
        manifest.append(dict(case=name,npz_sha256=sha(folder/f'{name}.npz'),json_sha256=sha(folder/f'{name}.json')));print('demo',name,flush=True)
    dump(folder/'manifest.json',dict(cases=manifest,wall_s=time.perf_counter()-start,freeze_sha256=sha(OUT/'results/model_freeze.json'),protocol_sha256=sha(ROOT/'configs/demo_protocol.json')))
if __name__=='__main__':main()

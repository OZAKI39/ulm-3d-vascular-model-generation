"""Validation selection, separate group calibration, then locked final evaluation."""
import sys,time,json,argparse,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,torch
from torch.utils.data import DataLoader
from inference import Predictor
from learning import CachedTriples
from vcommon import ROOT,OUT,DATA,P,LIMITS,NAMES,dump,sha,identity,stats,clean

@torch.inference_mode()
def predict(checkpoint,split,corruption=None,batch=128):
    model=Predictor(checkpoint);dl=DataLoader(CachedTriples(split,corruption),batch_size=batch,num_workers=4,pin_memory=True)
    pieces=[];start=time.perf_counter()
    for b in dl:
        raw=model.batch(b['images'],b['dt']);out=model.outputs(raw,b['dt'].numpy());pieces.append(out)
    arrays={k:np.concatenate([r[k] for r in pieces]) for k in pieces[0]}
    return arrays,dict(wall_s=time.perf_counter()-start,weight_sha256=model.sha,parameters=model.parameter_count,name=model.name)

def choose_threshold(score,bad,eligible,minimum=20):
    order=np.argsort(score);order=order[eligible[order]&np.isfinite(score[order])]
    if not len(order):return None,dict(n=0,failed=0,reason='NO_ELIGIBLE')
    # Evaluate only complete tied-score prefixes to preserve executable threshold semantics.
    risk=np.cumsum(bad[order])/np.arange(1,len(order)+1);ends=np.r_[score[order][1:]!=score[order][:-1],True]
    feasible=np.where((risk<=P['selection']['target_validation_accepted_failure_fraction'])&(np.arange(len(order))+1>=minimum)&ends)[0]
    if not len(feasible):return None,dict(n=0,failed=0,reason='NO_FEASIBLE_VALIDATION_RISK_THRESHOLD')
    k=feasible[-1];return float(score[order[k]]),dict(n=int(k+1),failed=int(bad[order[:k+1]].sum()),risk=float(risk[k]))

def freeze():
    if (OUT/'results/model_freeze.json').exists():raise RuntimeError('MODEL_FREEZE_ALREADY_EXISTS')
    records={};champions={}
    for folder in sorted((OUT/'training').glob('*_seed*')):
        if not (folder/'summary.json').exists():continue
        checkpoint=folder/'best.pt';name=folder.name
        val,timing=predict(checkpoint,'validation');truth=np.load(DATA/'validation/truth.npy');label=np.load(DATA/'validation/label.npy')
        err=np.abs(val['mean']-truth);bad=(label!=0)|np.any(err>LIMITS,axis=1)|~val['algebra_valid']
        eligible=(val['probability'].argmax(1)==0)&val['algebra_valid']
        threshold,selection=choose_threshold(val['score'],bad,eligible)
        field_thresholds=[];field_selection=[]
        for j in range(6):
            th,sel=choose_threshold(val['score'],(label!=0)|(err[:,j]>LIMITS[j])|~val['algebra_valid'],eligible);field_thresholds.append(th);field_selection.append(sel)
        cal,ctime=predict(checkpoint,'calibration');cy=np.load(DATA/'calibration/truth.npy');cl=np.load(DATA/'calibration/label.npy')
        metadata=[json.loads(x) for x in (DATA/'calibration/index.jsonl').read_text().splitlines()]
        groups=sorted(set((x['source_frame'],x['case']) for x in metadata if x['class_id']==0))
        errors=np.abs(cal['mean']-cy)/np.maximum(cal['sigma'],1e-9);scores=[]
        for frame,case in groups:
            ids=[i for i,m in enumerate(metadata) if m['source_frame']==frame and m['case']==case and m['class_id']==0]
            scores.append(np.max(errors[ids],axis=0))
        scores=np.array(scores);quantiles={}
        for level in P['calibration']['nominal_levels']:
            k=math.ceil((len(groups)+1)*level)
            if k>len(groups):quantiles[str(level)]=[None]*6
            else:quantiles[str(level)]=np.sort(scores,axis=0)[k-1].tolist()
        record=dict(checkpoint=str(checkpoint.relative_to(OUT)),weight_sha256=sha(checkpoint),threshold=threshold,selection=selection,field_thresholds=field_thresholds,field_selection=field_selection,quantiles=quantiles,calibration_groups=len(groups),group_definition='source frame x motion case (one recording; groups not independent experiments)',validation_timing=timing,calibration_timing=ctime,validation_score=json.loads((folder/'summary.json').read_text())['best_validation'])
        records[name]=record
        np.savez_compressed(OUT/'results'/f'{name}_validation.npz',**val)
        np.savez_compressed(OUT/'results'/f'{name}_calibration.npz',**cal)
        base=json.loads((folder/'config.json').read_text())['name']
        if base not in champions or record['validation_score']<records[champions[base]]['validation_score']:champions[base]=name
        dump(OUT/'results/calibration_in_progress.json',records);print('selected/calibrated',name,selection,flush=True)
    dump(OUT/'results/model_freeze.json',dict(identity=identity(),models=records,champions_by_validation_only=champions,data_manifest_hash={s:sha(DATA/s/'manifest.json') for s in ['train','validation','calibration','final_test']},final_test_evaluated=False))

def evaluate():
    frozen=json.loads((OUT/'results/model_freeze.json').read_text())
    for split,hashval in frozen['data_manifest_hash'].items():assert sha(DATA/split/'manifest.json')==hashval
    if (OUT/'results/final_evaluation_complete.json').exists():raise RuntimeError('FINAL_TEST_ALREADY_EVALUATED; do not tune/retest silently')
    rows=[]
    for name,r in frozen['models'].items():
        cp=OUT/r['checkpoint'];assert sha(cp)==r['weight_sha256']
        pred,timing=predict(cp,'final_test');file=OUT/'results'/f'{name}_final_test.npz';np.savez_compressed(file,**pred);rows.append(dict(model=name,corruption=None,file=file.name,sha256=sha(file),timing=timing))
        if name in [frozen['champions_by_validation_only'].get(k) for k in ['direct_fused','direct_fused_robust','direct_heat']]:
            for corruption in ['drop','false','jitter','severe']:
                pred,timing=predict(cp,'final_test',corruption);file=OUT/'results'/f'{name}_final_test_{corruption}.npz';np.savez_compressed(file,**pred);rows.append(dict(model=name,corruption=corruption,file=file.name,sha256=sha(file),timing=timing))
        print('final evaluated',name,flush=True)
    dump(OUT/'results/final_evaluation_complete.json',dict(freeze_sha256=sha(OUT/'results/model_freeze.json'),identity=identity(),evaluations=rows))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','test']);a=p.parse_args();torch.set_num_threads(4)
 (OUT/'results').mkdir(exist_ok=True)
 freeze() if a.stage=='freeze' else evaluate()

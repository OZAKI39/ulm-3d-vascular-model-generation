"""Development -> validation -> calibration -> immutable final evaluation."""
import sys,time,json,traceback,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from benchmark import *
import cv2,torch
MP=json.loads((ROOT/'configs/model_protocol.json').read_text())
METHODS=['PIV','CPD','POT_partial','POT_unbalanced','LRAFT_ST','CoTracker3','UnLiteFlowNet']
def load_split(split):
    index=json.loads((OUT/'data'/split/'index.json').read_text());det=old.CandidateDetector();result=[]
    for meta in index:
        path=OUT/meta['file'];assert sha(path)==meta['sha256'];d=dict(np.load(path));t=time.perf_counter();objects=[det(f)[0] for f in d['frames']]
        result.append((meta,d,objects,(time.perf_counter()-t)*1000))
    return result
def evaluate(split,name,model,config,samples,key):
    folder=OUT/'results/benchmark'/split/key;folder.mkdir(parents=True,exist_ok=True);rows=[];start=time.perf_counter()
    for meta,d,objects,detect_ms in samples:
        path=folder/f'{meta["index"]:04d}.json'
        if path.exists():
            row=json.loads(path.read_text());assert row['input_sha256']==meta['sha256'] and row['config']==config;rows.append(row);continue
        try:r=predict(name,model,d['frames'],objects,d['times'],config)
        except Exception as exc:r=dict(estimate=old.empty('METHOD_EXCEPTION',name),quality_score=0.,timing={},failure=repr(exc),traceback=traceback.format_exc(),vectors=[])
        score=scored(r,d['truth'],int(d['label']))
        row=dict(method=name,variant=key,config=config,index=meta['index'],case=meta['case'],source_frame=meta['source_frame'],seed=meta['seed'],split=split,input_sha256=meta['sha256'],truth=d['truth'],class_id=int(d['label']),**score,estimate=r['estimate'],quality_score=r['quality_score'],timing=r['timing'],detect_16_frames_ms=detect_ms,candidate_counts=[len(x) for x in objects],vectors=r.get('vectors',[]),failure=r.get('failure'))
        # Known marker-motion EPE uses generated moving locations, not background identity.
        if 'field' in r:
            from motion_models import sample_field
            xy=d['points'][-3];dd=d['points'][-1]-xy;pd=sample_field(r['field'],xy);good=np.isfinite(pd).all(1)&(xy[:,0]>=0)&(xy[:,0]<256)&(xy[:,1]>=0)&(xy[:,1]<256)
            row['moving_marker_epe_px']=float(np.linalg.norm(pd[good]-dd[good],axis=1).mean()) if good.any() else None
        dump(path,row);rows.append(clean(row))
    dump(folder/'COMPLETE.json',dict(identity=identity(),n=len(rows),config=config,wall_s=time.perf_counter()-start))
    return rows
def main():
    cv2.setNumThreads(1);torch.set_num_threads(2);start=time.perf_counter();monitor=ResourceMonitor();torch.cuda.reset_peak_memory_stats()
    freeze_path=OUT/'results/model_freeze.json';assert not freeze_path.exists(),'REFUSE_TO_RETUNE_FROZEN_MODEL'
    complete=json.loads((OUT/'data/DATASET_COMPLETE.json').read_text());dev=load_split('development');selected={};selection={}
    for name in METHODS:
        model=construct(name);variants=[{}]
        if name in ['LRAFT_ST','UnLiteFlowNet']:variants=[{'mode':m} for m in MP['model_variants']['dense_preprocess']]
        if name.startswith('POT_'):variants=[{'distance_scale':v} for v in MP['model_variants']['OT_distance_scales']]
        candidates=[]
        for i,config in enumerate(variants):
            rows=evaluate('development',name,model,config,dev,name+f'_v{i}')
            vals=[min(20,r['joint_normalized_error']) if r['joint_normalized_error'] is not None else 20 for r in rows if r['class_id']==0]
            candidates.append(dict(config=config,score_mean_clipped_joint_error=float(np.mean(vals)),variant=name+f'_v{i}'))
        winner=min(candidates,key=lambda v:v['score_mean_clipped_joint_error']);selected[name]=winner['config'];selection[name]=candidates
        print('development',name,winner,flush=True);del model;torch.cuda.empty_cache()
    dump(OUT/'results/development_selection.json',dict(identity=identity(),selected=selected,candidates=selection,rule='min mean clipped joint normalized raw error on global development examples; missing=20; no final/validation access',dataset_sha256=sha(OUT/'data/DATASET_COMPLETE.json')))
    del dev
    freeze=dict(identity=identity(),data=complete,selected=selected,limits=LIMITS,gates={},calibration={},scope=MP['scope'],true_real_accuracy='NOT_VERIFIED',weights={})
    for split in ['validation','calibration']:
        samples=load_split(split)
        for name in METHODS:
            model=construct(name);rows=evaluate(split,name,model,selected[name],samples,name)
            if model is not None:freeze['weights'][name]=dict(path=str(model.path),sha256=sha(model.path))
            if split=='validation':
                thresholds=sorted(set([0.]+[float(r['quality_score']) for r in rows]));options=[]
                for q in thresholds:
                    accepted=[r for r in rows if r['base_valid'] and r['quality_score']>=q];bad=sum(r['bad_if_accepted'] for r in accepted)
                    if len(accepted)>=10 and bad/len(accepted)<=.05:options.append(dict(threshold=q,n=len(accepted),bad=bad,risk=bad/len(accepted)))
                best=max(options,key=lambda a:a['n']) if options else dict(threshold=None,n=0,bad=0,risk=None)
                freeze['gates'][name]=dict(**best,rule='validation empirical joint risk <=.05, >=10 acceptances; no statistical guarantee')
                print('validation',name,freeze['gates'][name],flush=True)
            else:
                err=[r['absolute_error'] for r in rows if r['class_id']==0 and r['absolute_error'] is not None]
                freeze['calibration'][name]=dict(n=len(err),absolute_error_q90=np.percentile(err,90,axis=0) if err else None,scope='empirical synthetic residuals; not real confidence intervals; source replicates dependent')
            del model;torch.cuda.empty_cache()
        del samples
    dump(freeze_path,freeze);dump(OUT/'results/FINAL_OPEN_RECORD.json',dict(identity=identity(),freeze_sha256=sha(freeze_path),timestamp=time.time(),policy='final opened after frozen dev selection / validation risk gates / calibration'))
    samples=load_split('final_test')
    for name in METHODS:
        model=construct(name);rows=evaluate('final_test',name,model,selected[name],samples,name)
        threshold=freeze['gates'][name]['threshold'];out=[]
        for r in rows:
            r['risk_gate_valid']=r['base_valid'] and threshold is not None and r['quality_score']>=threshold
            r['published_estimate']=apply_gate(r['estimate'],r['quality_score'],threshold);out.append(r)
        dump(OUT/'results'/f'final_{name}.json',out);print('final',name,len(out),'base',sum(r['base_valid'] for r in out),'gate',sum(r['risk_gate_valid'] for r in out),flush=True);del model;torch.cuda.empty_cache()
    dump(OUT/'results/MODEL_BENCHMARK_COMPLETE.json',dict(identity=identity(),wall_s=time.perf_counter()-start,methods=METHODS,final_count=len(samples),resource=monitor.close(),peak_cuda_allocated_mb=torch.cuda.max_memory_allocated()/2**20,freeze_sha256=sha(freeze_path),training_executed=False))
if __name__=='__main__':main()

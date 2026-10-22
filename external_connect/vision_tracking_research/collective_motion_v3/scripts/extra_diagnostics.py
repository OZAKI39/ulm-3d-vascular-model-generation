"""Bounded real clip, historical diagnostics, official sample, solver warnings."""
import sys,time,json,contextlib,io,warnings,traceback,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from benchmark import *
import cv2,torch
DP=json.loads((ROOT/'configs/diagnostic_protocol.json').read_text())
def real_clip():
    freeze=json.loads((OUT/'results/model_freeze.json').read_text());begin,end=DP['selected_real_clip'];sources=old.read_frames(range(DP['history_begin'],end),VIDEO);det=old.CandidateDetector();images=[];objs=[];ts=[];x,y,w,h=P['roi_xywh'];start=time.perf_counter()
    for fid,(bgr,pts) in sources.items():o,g=det(bgr[y:y+h,x:x+w]);images.append(g);objs.append(o);ts.append(pts)
    allrows={}
    for name in DP['methods']:
        model=construct(name);rows=[]
        for fid in range(begin,end):
            j=fid-DP['history_begin'];r=predict(name,model,images[j-15:j+1],objs[j-15:j+1],ts[j-15:j+1],freeze['selected'][name]);e=r['estimate'];published=apply_gate(e,r['quality_score'],freeze['gates'][name]['threshold']) if name!='PIV' else e
            rows.append(dict(frame_id=fid,pts_s=ts[j],estimate=e,published_estimate=published,vectors=r['vectors'],quality_score=r['quality_score'],timing=r['timing'],objects=objs[j],context_frame_ids=list(range(fid-15,fid+1)) if name=='CoTracker3' else [fid-2,fid-1,fid]))
        allrows[name]=rows;del model;torch.cuda.empty_cache();print('real clip',name,len(rows),flush=True)
    dump(OUT/'results/real_comparison_clip.json',dict(identity=identity(),protocol=DP,source_sha256=sha(VIDEO),freeze_sha256=sha(OUT/'results/model_freeze.json'),methods=allrows,wall_s=time.perf_counter()-start,real_accuracy='NOT_VERIFIED',causal=True))
def historical():
    freeze=json.loads((OUT/'results/model_freeze.json').read_text());det=old.CandidateDetector();source=old.read_frames([840],VIDEO)[840];cases=json.loads((ROOT.parent/'collective_motion_v1/configs/protocol.json').read_text())['cases'];samples=[]
    for case in cases:
        im,meta=old.original_synthetic.build_sequence(case,7201,source,length=16);ts=np.arange(16)*meta['dt'];rel=meta['H'][15]@np.linalg.inv(meta['H'][13]);truth=old.parameters(rel[:2,:2],rel[:2,2],ts[-1]-ts[-3],np.array([128,128]));samples.append(('F10U_historical',case,im,[det(f)[0] for f in im],ts,truth,1 if case=='local_opposite' else 0))
    v=Path('/workspace/f10v_vision_20261017/data/final_test');index=[json.loads(l) for l in (v/'index.jsonl').read_text().splitlines()];raw=np.load(v/'raw.npy',mmap_mode='r');truths=np.load(v/'truth.npy');labels=np.load(v/'label.npy');dt=np.load(v/'dt.npy')
    for i in DP['historical_f10v_indices']:
        im=np.asarray(raw[i]);samples.append(('F10V_exposed_final_diagnostic',index[i]['case'],im,[det(f)[0] for f in im],np.r_[0,np.cumsum(dt[i])],truths[i],int(labels[i])))
    rows=[];start=time.perf_counter()
    for name in ['PIV','CPD','POT_partial','POT_unbalanced','LRAFT_ST','CoTracker3','UnLiteFlowNet']:
        model=construct(name)
        for scope,case,im,obj,ts,truth,label in samples:
            try:r=predict(name,model,im,obj,ts,freeze['selected'][name]);s=scored(r,truth,label);rows.append(dict(method=name,scope=scope,case=case,context_frames=len(im),truth=truth,class_id=label,**s,estimate=r['estimate'],quality_score=r['quality_score'],timing=r['timing']))
            except Exception as e:rows.append(dict(method=name,scope=scope,case=case,status='FAILED',error=repr(e)))
        del model;torch.cuda.empty_cache();print('historical',name,flush=True)
    dump(OUT/'results/historical_diagnostics.json',dict(identity=identity(),rows=rows,wall_s=time.perf_counter()-start,scope='EXPOSED_HISTORICAL_DIAGNOSTICS_ONLY',f10v_raw_sha256=sha(v/'raw.npy'),f10v_index_sha256=sha(v/'index.jsonl'),no_retuning=True))
def ot_audit():
    from motion_models import optimal_transport
    freeze=json.loads((OUT/'results/model_freeze.json').read_text());det=old.CandidateDetector();index=json.loads((OUT/'data/final_test/index.json').read_text());rows=[]
    for meta in index:
        d=dict(np.load(OUT/meta['file']));objects=[det(f)[0] for f in d['frames']]
        for method in ['POT_partial','POT_unbalanced']:
            buf=io.StringIO()
            with contextlib.redirect_stdout(buf),warnings.catch_warnings(record=True) as warn:
                warnings.simplefilter('always');r=optimal_transport(d['frames'],objects,d['times'],method=method[4:],**freeze['selected'][method])
            messages=r.get('solver_warning','')+buf.getvalue()+'\n'+'\n'.join(str(w.message) for w in warn);prior=json.loads((OUT/'results'/f'final_{method}.json').read_text())[meta['index']];rows.append(dict(index=meta['index'],method=method,case=meta['case'],base_valid=all(r['estimate']['valid'].values()),previous_base_valid=prior['base_valid'],warning=messages.strip(),transported_mass=r.get('transported_mass'),solver_log=r.get('solver_log'),estimate=r['estimate']))
    dump(OUT/'results/ot_numerical_audit.json',dict(identity=identity(),rows=rows,numerical_warning_and_accepted=[dict(index=r['index'],method=r['method']) for r in rows if r['warning'] and r['base_valid']],scope='Exact frozen final inputs and settings, numerical audit only; no retuning'))
def official_example():
    from PIL import Image
    import flowiz
    model=construct('UnLiteFlowNet');folder=THIRD/'UnLiteFlowNet-PIV/sample_data';rows=[];start=time.perf_counter()
    for p in sorted(folder.glob('*/*img1.tif')):
        q=p.with_name(p.name.replace('img1','img2'));gt=p.with_name(p.name.replace('img1.tif','flow.flo'));a=np.asarray(Image.open(p));b=np.asarray(Image.open(q));truth=flowiz.read_flow(str(gt));f,timing=model.flow(a,b,'raw');rows.append(dict(sample=str(p.relative_to(folder)),input_sha256=[sha(p),sha(q),sha(gt)],epe_px=float(np.linalg.norm(f-truth,axis=2).mean()),timing=timing))
    dump(OUT/'results/UnLite_official_examples.json',dict(identity=identity(),rows=rows,wall_s=time.perf_counter()-start,weight_sha256=sha(model.path),scope='Official bundled sample_data, official estimate() and flowiz GT reader; deterministic all examples instead of main.py random plotting selector'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('task',choices=['real','historical','ot','official']);a=p.parse_args();cv2.setNumThreads(1);torch.set_num_threads(2)
    {'real':real_clip,'historical':historical,'ot':ot_audit,'official':official_example}[a.task]()

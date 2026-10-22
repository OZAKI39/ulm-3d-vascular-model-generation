"""Bounded bulk proposals from published Cellpose/µSAM models, never ground truth."""
from pathlib import Path
import os,sys,json,time,traceback,select,hashlib,subprocess
import importlib.metadata
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MICROSAM_CACHEDIR',str(ROOT/'models/microsam_cache'))
os.environ.setdefault('CELLPOSE_LOCAL_MODELS_PATH',str(ROOT/'models/cellpose4'))
os.environ.setdefault('OMP_NUM_THREADS','2')
import numpy as np
from PIL import Image
from annotation_store import sha
from annotation_regions import normalize_box
from annotation_sizes import size_plan,collect_candidates
from sam_auto_training import gpu_slot

MODEL=None;MODEL_NAME=None;MODEL_SHA=None

def emit(**message):print(json.dumps(message,ensure_ascii=False),flush=True)

def get_model(name):
    global MODEL,MODEL_NAME,MODEL_SHA
    if name==MODEL_NAME:return MODEL
    import torch,gc
    torch.set_num_threads(2)
    if not torch.cuda.is_available():raise RuntimeError('GPU unavailable; existing labels are unchanged')
    MODEL=None;MODEL_NAME=None;MODEL_SHA=None;gc.collect();torch.cuda.empty_cache()
    if name=='microsam_lm':
        from micro_sam.instance_segmentation import get_predictor_and_decoder,InstanceSegmentationWithDecoder
        predictor,decoder=get_predictor_and_decoder('vit_b_lm',device='cuda')
        MODEL=InstanceSegmentationWithDecoder(predictor,decoder)
        paths=[ROOT/'models/microsam_cache/models'/p for p in ('vit_b_lm','vit_b_lm_decoder')]
    elif name in ('cellpose4','cyto3'):
        from cellpose import models
        if name=='cellpose4':MODEL=models.CellposeModel(gpu=True,pretrained_model='cpsam_v2')
        else:
            checkpoint=ROOT/'models/cellpose/cyto3'
            if not checkpoint.is_file():raise FileNotFoundError('Missing required cyto3 checkpoint: '+str(checkpoint))
            MODEL=models.CellposeModel(gpu=True,pretrained_model=str(checkpoint))
        paths=[Path(MODEL.pretrained_model)]
    else:raise ValueError('Unknown proposal model')
    MODEL_SHA={str(p.relative_to(ROOT)):sha(p) for p in paths};MODEL_NAME=name
    return MODEL

def propose(req):
    start=time.monotonic();image_path=Path(req['image'])
    if sha(image_path)!=req['image_sha256']:raise ValueError('Image SHA mismatch')
    with Image.open(image_path) as im:image=np.asarray(im.convert('RGB'))
    roi=normalize_box(req['roi_yxyx'],image.shape[:2])
    y0,x0,y1,x1=roi
    if max(y1-y0,x1-x0)>768 or min(y1-y0,x1-x0)<32:raise ValueError('请选择边长 32–768 px 的小区域')
    plan=size_plan(req)
    pad=32;cy0=max(0,y0-pad);cx0=max(0,x0-pad);cy1=min(image.shape[0],y1+pad);cx1=min(image.shape[1],x1+pad)
    crop=image[cy0:cy1,cx0:cx1]
    if req.get('invert',False):crop=255-crop
    emit(event='progress',token=req['token'],status='加载模型 / 等待 GPU，已有标注保留')
    with gpu_slot(interactive=True):
        import torch
        torch.cuda.reset_peak_memory_stats()
        model=get_model(req['model']);loaded=time.monotonic()
        raw_by_scale=[];scale_timings=[]
        for index,d in enumerate(plan['inference_diameters_px']):
            stage=time.monotonic();detail=f'参考检测直径 {d:.1f} px' if d is not None else '模型原生尺度'
            emit(event='progress',token=req['token'],status=f'检测 {index+1}/{len(plan["inference_diameters_px"])}：{detail}')
            if req['model']=='microsam_lm':
                model.initialize(crop);raw=model.generate(min_size=5)
            else:
                kw={'batch_size':1,'min_size':5,'flow_threshold':.4,'cellprob_threshold':0,'augment':False}
                if req['model']=='cyto3':kw.update(channels=[0,0],diameter=d)
                raw,_,_=model.eval(crop,**kw)
            torch.cuda.synchronize()
            raw_by_scale.append(np.asarray(raw,dtype=np.uint32))
            scale_timings.append({'scale_index':index,'diameter_px':d,'inference_s':time.monotonic()-stage})
        torch.cuda.synchronize();peak=torch.cuda.max_memory_allocated()/2**20
    labels,objects,rejected=collect_candidates(raw_by_scale,[y0-cy0,x0-cx0,y1-cy0,x1-cx0],roi,plan)
    out=Path(req['output']);out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,labels=labels,raw_labels=raw_by_scale[0],
                        **{f'raw_labels_s{i:02d}':a for i,a in enumerate(raw_by_scale)})
    execution={'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
               'python':sys.executable,'source_sha256':{p:sha(ROOT/p) for p in ['scripts/preannotation_worker.py','src/annotation_regions.py','src/annotation_sizes.py']},
               'model_package_version':importlib.metadata.version('micro-sam' if req['model']=='microsam_lm' else 'cellpose')}
    result={'token':req['token'],'image_sha256':req['image_sha256'],'model':req['model'],'model_sha256':MODEL_SHA,'execution':execution,
            'roi_yxyx':roi,'context_yxyx':[cy0,cx0,cy1,cx1],'objects':objects,'rejected':rejected,
            'status':'UNREVIEWED_PROPOSALS','load_s':loaded-start,'inference_and_filter_s':time.monotonic()-loaded,
            'worker_wall_s':time.monotonic()-start,'peak_cuda_mib':peak,'output_sha256':sha(out),
            'size_plan':plan,'scale_timings':scale_timings,
            'merge_policy':'whole_mask_scale_match_ranking_no_clipping',
            'input_inverted':bool(req.get('invert',False)),'automatic_acceptance':False,'weights_are_experiment_calibrated':False}
    if plan['mode']=='DIAMETER_RANGE':result['diameter_range_px']=plan['range_px']
    else:result['diameter_px']=plan['nominal_diameter_px']
    out.with_suffix('.json').write_text(json.dumps(result,indent=2))
    return result

def main():
    if len(sys.argv)>1 and sys.argv[1]!='--server':
        req=json.loads(Path(sys.argv[1]).read_text());emit(event='result',**propose(req));return
    while select.select([sys.stdin],[],[],120)[0]:
        line=sys.stdin.readline()
        if not line:return
        req={}
        try:req=json.loads(line);emit(event='result',**propose(req))
        except Exception:
            error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
            emit(event='error',token=req.get('token'),error=error)

if __name__=='__main__':main()

"""Finite Cellpose/OpenCV adaptation on a frozen new-microscope annotation set.

Test images are first decoded after validation selection is saved. The Cellpose
proposal path uses exactly the annotation worker's padding, scales and merge rule.
"""
from pathlib import Path
import argparse, collections, hashlib, importlib.metadata, json, logging
import shutil, signal, subprocess, sys, time, traceback

BASE = Path(__file__).resolve().parents[1]
LIVE = BASE.parent/'realtime_collective_observation_v1'
sys.path.insert(0, str(BASE/'src'))
sys.path.insert(0, str(LIVE/'src'))
import cv2
import numpy as np
from annotation_store import CLASSES, sha
from annotation_snapshot import freeze_reviewed, export_frozen
from annotation_regions import reviewed_boxes
from annotation_sizes import size_plan, collect_candidates
from vision_evaluation import evaluate_instances, aggregate
from microscope_detector import LEGACY, detect


def dump(p, value):
    Path(p).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False))


def run(root, out):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    cv2.setNumThreads(2)
    execution = dict(head=subprocess.check_output(['git','rev-parse','HEAD'], cwd=BASE, text=True).strip(),
                     sources={str(p):sha(p) for p in [Path(__file__), BASE/'src/vision_evaluation.py',
                     BASE/'src/annotation_snapshot.py', BASE/'src/annotation_sizes.py', LIVE/'src/microscope_detector.py']},
                     python=sys.executable, cellpose_version=importlib.metadata.version('cellpose'))
    source_dir=out/'source_snapshot';source_dir.mkdir()
    for source in execution['sources']:shutil.copy2(source,source_dir/Path(source).name)
    freeze = out/'freeze'; frozen = freeze_reviewed(root, freeze, execution); export_frozen(freeze)
    summary = collections.defaultdict(lambda:dict(frames=0, single_bubbles=0))
    for row in frozen['rows']:
        s=summary[row['record']['split']];s['frames']+=1
        s['single_bubbles']+=row['annotation']['reviewed_instance_counts']['SINGLE_BUBBLE']
    dump(out/'dataset.json', dict(splits=dict(summary), exclusions=frozen['excluded']))
    train_export=json.loads((freeze/'cellpose_train/export_manifest.json').read_text())
    import tifffile
    train_images=[];train_labels=[];diameters=[]
    for item in train_export:
        if item['status']!='EXPORTED':continue
        folder=freeze/'cellpose_train'
        train_images.append(cv2.imread(str(folder/(item['id']+'.png'))))
        labels=tifffile.imread(folder/(item['id']+'_masks.tif'));train_labels.append(labels)
        diameters += [float(2*np.sqrt(np.count_nonzero(labels==n)/np.pi)) for n in np.unique(labels) if n]
    if not train_images:raise ValueError('No eligible reviewed training crops')
    # Fixed, small grid, includes the exact current LIVE detector as a baseline.
    grid=[dict(LEGACY)]
    for sigma in (2.,4.,6.):
        for threshold in (4.,8.,12.):
            for area in (5,8):
                candidate=dict(LEGACY,sigma_px=sigma,threshold=threshold,min_area_px=area)
                if candidate not in grid:grid.append(candidate)
    protocol=dict(scope='REVIEWED_ROIS_SAME_RECORDING_ONLY', execution=execution,
        freeze_sha256=sha(freeze/'frozen.json'),train_label_median_diameter_px=float(np.median(diameters)),
        evaluation=dict(center_radius_px=6,border_margin_px=6,mask_iou=.5,
            large_predictions='Count as predictions; no silent area exclusion',
            boundary_truth='Exclude truncated or near-edge truth and overlapping predictions; report counts'),
        opencv_grid=grid,opencv_selection='Best training center F1, then fewer FP, then legacy preference',
        validation_admission='Center F1 improvement >=0.02, recall >=baseline, precision >=baseline-0.03',
        cellpose=dict(model='cyto3',epochs=100,learning_rate=.1,SGD=True,batch_size=2,bsize=224,
            weight_decay=1e-4,max_training_wall_s=300,seed=20261019,validation_to_optimizer=False,
            checkpoint_selection='Final epoch only',scales_px=[8,13.856406460551018,24],
            context_pad_px=32,diameter_filter='REFERENCE_ONLY',min_size=5),
        test_gate='Chosen variant F1 >= baseline-0.02 and precision >= baseline-0.03; no reselection on test',
        max_total_wall_s=900,channel_mask='NOT_DRAWN',live_camera_test='NOT_RUN',
        test_usage='Only after selection.json is frozen; no hyperparameter changes after test',
        prohibited_claims=['cross-session validation','physical bubble identity','PIV velocity accuracy'])
    dump(out/'protocol.json',protocol)
    all_rows=[];frames_cache={};pred_dir=out/'predictions';pred_dir.mkdir()
    def frames(split):
        if split=='test' and not (out/'selection.json').exists():raise RuntimeError('Test locked until selection')
        if split not in frames_cache:
            frames_cache[split]=[]
            for row in frozen['rows']:
                if row['record']['split']!=split:continue
                with np.load(freeze/row['frozen_mask_file']) as z:masks={k:z[k].copy() for k in CLASSES}
                frames_cache[split].append(dict(row=row,masks=masks,image=cv2.imread(str(freeze/row['frozen_image_file'])),
                                                boxes=reviewed_boxes(row['annotation'],masks)))
        return frames_cache[split]
    def evaluate(fn,name,split,persist=True):
        rows=[]
        for f in frames(split):
            full=f['masks']['SINGLE_BUBBLE'];counts=dict(zip(*np.unique(full,return_counts=True)))
            for i,box in enumerate(f['boxes']):
                if time.monotonic()-start>900:raise TimeoutError('Total study budget exceeded')
                tick=time.monotonic();labels,objects=fn(f['image'],box);ms=(time.monotonic()-tick)*1000
                y0,x0,y1,x1=box
                metrics=evaluate_instances(full[y0:y1,x0:x1],labels,objects,
                    (f['masks']['UNCERTAIN'][y0:y1,x0:x1]>0)|(f['masks']['OUT_OF_FOCUS'][y0:y1,x0:x1]>0),
                    counts,max_prediction_area=None)
                key=f'{name}_{f["row"]["record"]["id"]}_roi{i}'
                row=dict(model=name,split=split,image_id=f['row']['record']['id'],roi_index=i,
                         box_yxyx=box,detect_ms=ms,metrics=metrics,artifact=key)
                if persist:np.savez_compressed(pred_dir/(key+'.npz'),labels=labels,truth=full[y0:y1,x0:x1])
                rows.append(row)
        all_rows.extend(rows);return aggregate(rows)
    def ocv(config):
        def fn(image,box):
            y0,x0,y1,x1=box;_,labels,objects=detect(cv2.cvtColor(image[y0:y1,x0:x1],cv2.COLOR_BGR2GRAY),config)
            # Only accepted components remain in the displayed/evaluated map.
            keep=np.zeros(int(labels.max())+1,bool)
            for o in objects:keep[o['label']]=True
            return np.where(keep[labels],labels,0),objects
        return fn
    scores=[evaluate(ocv(c),f'opencv_grid_{i}','train',False) for i,c in enumerate(grid)]
    best=max(range(len(grid)),key=lambda i:(scores[i]['center']['f1'] or 0,-scores[i]['center']['fp'],-i))
    dump(out/'opencv_search.json',dict(configs=grid,scores=scores,selected_index=best,selection_split='train'))
    metrics={}
    for name,c in [('opencv_original',LEGACY),('opencv_calibrated',grid[best])]:
        for split in ('train','validation'):metrics[f'{name}/{split}']=evaluate(ocv(c),name,split)
    print('OpenCV training search and validation complete',flush=True)
    import torch
    from cellpose import models,train
    torch.set_num_threads(2);torch.manual_seed(20261019);np.random.seed(20261019)
    base_weights=(root/'models/cellpose/cyto3').resolve()
    model=models.CellposeModel(gpu=True,pretrained_model=str(base_weights))
    plan=size_plan(dict(model='cyto3',diameter_range_px=[8,24],diameter_filter='REFERENCE_ONLY'))
    def cp(image,box):
        y0,x0,y1,x1=box;cy0=max(0,y0-32);cx0=max(0,x0-32);cy1=min(image.shape[0],y1+32);cx1=min(image.shape[1],x1+32)
        crop=cv2.cvtColor(image[cy0:cy1,cx0:cx1],cv2.COLOR_BGR2RGB);raw=[]
        for d in plan['inference_diameters_px']:
            labels,_,_=model.eval(crop,channels=[0,0],diameter=d,batch_size=1,min_size=5,
                                  flow_threshold=.4,cellprob_threshold=0,augment=False)
            raw.append(np.asarray(labels,np.uint32))
        torch.cuda.synchronize()
        labels,objs,_=collect_candidates(raw,[y0-cy0,x0-cx0,y1-cy0,x1-cx0],box,plan)
        return labels,[dict(label=o['id'],area=o['area_px'],xy_local=[o['centroid_yx'][1]-x0,o['centroid_yx'][0]-y0]) for o in objs]
    # Original and tuned paths match the current multi-scale proposal interface.
    for split in ('train','validation'):metrics[f'cellpose_original/{split}']=evaluate(cp,'cellpose_original',split)
    dump(out/'pretraining_metrics.json',metrics)
    tick=time.monotonic()
    def timeout(*_):raise TimeoutError('Training budget 300 seconds exceeded')
    prev=signal.signal(signal.SIGALRM,timeout);signal.alarm(300)
    try:
        weights,losses,_=train.train_seg(model.net,train_data=train_images,train_labels=train_labels,
            channels=[0,0],normalize=True,n_epochs=100,batch_size=2,learning_rate=.1,SGD=True,
            weight_decay=1e-4,bsize=224,min_train_masks=5,save_path=out,save_every=100,
            model_name='cyto3_new_microscope_100ep')
        training=dict(status='COMPLETED',wall_s=time.monotonic()-tick,weights=str(weights),
                      weights_sha256=sha(weights),losses=losses.tolist(),train_crops=len(train_images))
    except Exception:
        training=dict(status='FAILED',wall_s=time.monotonic()-tick,error=traceback.format_exc())
    finally:signal.alarm(0);signal.signal(signal.SIGALRM,prev)
    dump(out/'training.json',training);print('Cellpose training '+training['status'],flush=True)
    if training['status']=='COMPLETED':
        model=models.CellposeModel(gpu=True,pretrained_model=training['weights'])
        for split in ('train','validation'):metrics[f'cellpose_tuned/{split}']=evaluate(cp,'cellpose_tuned',split)
    def improves(candidate,baseline):
        c=metrics[candidate+'/validation']['center'];b=metrics[baseline+'/validation']['center']
        return bool(c['f1']>=b['f1']+.02 and c['recall']>=b['recall'] and c['precision']>=b['precision']-.03)
    selected_ocv='opencv_calibrated' if improves('opencv_calibrated','opencv_original') else 'opencv_original'
    selected_cp='cellpose_tuned' if training['status']=='COMPLETED' and improves('cellpose_tuned','cellpose_original') else 'cellpose_original'
    selection=dict(opencv=selected_ocv,opencv_parameters=grid[best] if selected_ocv.endswith('calibrated') else LEGACY,
                   cellpose=selected_cp,cellpose_weights=training['weights'] if selected_cp.endswith('tuned') else str(base_weights),
                   selected_on='validation',test_accessed=False,selected_at_unix_s=time.time())
    dump(out/'selection.json',selection);selection_sha=sha(out/'selection.json')
    print('Selection frozen; now evaluating held-out test regions',flush=True)
    for name,params in [('opencv_original',LEGACY),(selected_ocv,selection['opencv_parameters'])]:
        if name+'/test' not in metrics:metrics[name+'/test']=evaluate(ocv(params),name,'test')
    for name,weights in [('cellpose_original',str(base_weights)),(selected_cp,selection['cellpose_weights'])]:
        if name+'/test' not in metrics:
            model=models.CellposeModel(gpu=True,pretrained_model=weights)
            metrics[name+'/test']=evaluate(cp,name,'test')
    assert sha(out/'selection.json')==selection_sha
    def heldout_gate(name,baseline):
        c=metrics[name+'/test']['center'];b=metrics[baseline+'/test']['center']
        return bool(c['f1'] is not None and c['f1']>=b['f1']-.02 and c['precision']>=b['precision']-.03)
    # Full-image inference cost only; no unreviewed whole-image accuracy claim.
    image=frames('train')[0]['image'];gray=cv2.cvtColor(image,cv2.COLOR_BGR2GRAY);times=[]
    for _ in range(30):
        t=time.monotonic();detect(gray,selection['opencv_parameters']);times.append(1000*(time.monotonic()-t))
    bench=dict(opencv_full_image_ms=dict(p50=float(np.median(times)),p95=float(np.percentile(times,95))),
               image_wh=[image.shape[1],image.shape[0]],scope='Detector only; not end-to-end LIVE fps')
    result=dict(status='COMPLETED',metrics=metrics,selection=selection,selection_sha256=selection_sha,
        opencv_update_admitted=selected_ocv=='opencv_calibrated' and heldout_gate(selected_ocv,'opencv_original'),
        cellpose_update_admitted=selected_cp=='cellpose_tuned' and heldout_gate(selected_cp,'cellpose_original'),
        benchmark=bench,total_wall_s=time.monotonic()-start,training_status=training['status'],
        scope='Local reviewed regions of one recording only',channel_mask='NOT_DRAWN',live_accuracy='NOT_MEASURED',
        training_ready=False,control_permission=False)
    dump(out/'rows.json',all_rows);dump(out/'summary.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='metrics'},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();logging.basicConfig(level=logging.INFO)
    try:run(a.workspace.resolve(),a.output.resolve())
    except Exception:
        if a.output.exists():dump(a.output/'FAILED.json',dict(error=traceback.format_exc()))
        raise

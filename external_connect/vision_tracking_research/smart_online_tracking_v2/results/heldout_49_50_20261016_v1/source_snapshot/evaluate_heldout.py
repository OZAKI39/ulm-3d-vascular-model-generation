"""Frozen-parameter test on already reviewed held-out ROIs. No training or tuning."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, json, sys, time, signal, traceback
import numpy as np
import cv2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from annotation_store import sha, CLASSES
from annotation_regions import reviewed_boxes, reviewed_counts
from vision_engine import VisionEngine
from vision_evaluation import evaluate_instances, aggregate
from runtime import identity
from sam_auto_training import snapshot as sam_training_snapshot

BASE = ROOT/'results/supervised_20261016_v1'
EXPECTED_CYTO3 = '2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c'


def write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False))


def evaluate(freeze, out):
    out.mkdir(exist_ok=False)
    started = time.monotonic()
    frozen = json.loads((freeze/'frozen.json').read_text())
    for name, digest in frozen['snapshot_sha256'].items():
        if sha(freeze/name) != digest:
            raise ValueError('Frozen input changed: '+name)
    selected_ids = ['video_000660', 'video_000678']
    selected = [r for r in frozen['rows'] if r['record']['id'] in selected_ids]
    if [r['record']['id'] for r in selected] != selected_ids:
        raise ValueError('Both requested test frames must have valid reviewed coverage')
    configs = {r: json.loads((BASE/f'selected_{r}.json').read_text()) for r in ['opencv','cellpose']}
    baseline = json.loads((BASE/'preregistered.json').read_text())
    if configs['opencv'] != {'contrast_threshold':10} or configs['cellpose'] != {
            'diameter':baseline['cellpose']['diameter_from_train_median_px'], 'min_size':5}:
        raise ValueError('Selected baseline changed; no test-driven model selection allowed')
    if sha(ROOT/'models/cellpose/cyto3') != EXPECTED_CYTO3:
        raise ValueError('Pretrained model checksum changed')
    for name in ['src/vision_engine.py','src/vision_evaluation.py']:
        if sha(ROOT/name) != baseline['execution']['source_sha256'][name]:
            raise ValueError('Evaluation/inference changed since validation: '+name)
    # Test masks are also barred from the existing background SAM training service.
    auto = sam_training_snapshot(ROOT)
    assert all(r['record']['split'] in ('train','validation') for r in auto['rows'])
    frames = []; coverage = []
    for row in selected:
        if row['record']['split'] != 'test' or row['annotation']['split'] != 'test':
            raise ValueError('Split identity mismatch')
        with np.load(freeze/row['frozen_mask_file']) as z:
            masks = {k:z[k].copy() for k in CLASSES}
        boxes = reviewed_boxes(row['annotation'],masks)
        # Reject overlapping review boxes instead of counting the same pixels twice.
        seen = np.zeros(masks['SINGLE_BUBBLE'].shape,bool)
        for y0,x0,y1,x1 in boxes:
            if seen[y0:y1,x0:x1].any():
                raise ValueError('Overlapping test ROIs require an explicit policy before scoring')
            seen[y0:y1,x0:x1] = True
        image = cv2.imread(str(freeze/row['frozen_image_file']))
        frames.append(dict(row=row,masks=masks,boxes=boxes,image=image))
        counts = reviewed_counts(row['annotation'],masks)
        coverage.append(dict(image_id=row['record']['id'],revision=row['annotation']['revision'],
                             review_scope=row['annotation']['review_scope'],boxes_yxyx=boxes,
                             reviewed_counts=counts,saved_counts=row['annotation']['instance_counts'],
                             unreviewed_saved_single_masks=row['annotation']['instance_counts']['SINGLE_BUBBLE']-counts['SINGLE_BUBBLE']))
    prereg = dict(utc=datetime.now(timezone.utc).isoformat(),source='SAVED_VIDEO_HELDOUT_TEST',
        selected_ids=selected_ids,freeze=str(freeze),freeze_sha256=sha(freeze/'frozen.json'),
        baseline_result=str(BASE),configs=configs,model_sha256=EXPECTED_CYTO3,
        config_sha256={r:sha(BASE/f'selected_{r}.json') for r in configs},coverage=coverage,
        primary='Same reviewed original-pixel ROI crops; center precision/recall/F1 at 6 px; mask F1 at IoU >= 0.5',
        secondary='Full-frame inference, scored only in the same reviewed ROIs; not full-frame accuracy',
        exclusion='Inherited 6 px center border and clipped truth; ignored predictions overlap excluded truth >=50%',
        no_test_training=True,no_test_tuning=True,no_test_model_selection=True,
        rejected_finetuned_model='NOT_RERUN; decision fixed by preceding validation',
        groups_and_tracks='NOT_EVALUATED_NO_LABELS',cross_session_generalization='NOT_VERIFIED',
        full_frame_human_review=False,wall_budget_s=180,execution=identity(ROOT),
        auto_sam_training_ids=[r['record']['id'] for r in auto['rows']],
        auto_sam_test_labels_excluded=True,
        figures='All 5 ROI examples and test/validation metrics; no outcome-based sample selection')
    write(out/'preregistered.json',prereg)
    src = out/'source_snapshot';src.mkdir()
    for name in ['scripts/evaluate_heldout.py','src/vision_engine.py','src/vision_evaluation.py','src/annotation_regions.py']:
        target=src/Path(name).name;target.write_bytes((ROOT/name).read_bytes())
    raw = out/'predictions';raw.mkdir()
    rows = []; timings = []
    oldfreeze = ROOT/'results/expert_freezes/vision2_supervised_20261016_v1'
    oldrows = json.loads((oldfreeze/'frozen.json').read_text())['rows']
    train = next(r for r in oldrows if r['record']['split']=='train')
    warm = cv2.imread(str(oldfreeze/train['frozen_image_file']))[322:491,679:882]
    for route in ['opencv','cellpose']:
        tick=time.monotonic();engine=VisionEngine(route,detector_config=configs[route])
        engine.detect(warm)
        timings.append(dict(route=route,scope='TRAIN_ONLY_WARMUP_AND_LOAD',wall_s=time.monotonic()-tick))
        for f in frames:
            image_id=f['row']['record']['id'];full=f['masks']['SINGLE_BUBBLE']
            counts={int(n):int(np.count_nonzero(full==n)) for n in np.unique(full) if n}
            for domain in ['ROI','FULL_FRAME']:
                if domain=='FULL_FRAME':
                    tick=time.monotonic();objects,_,labels=engine.detect(f['image'])
                    elapsed=time.monotonic()-tick
                    timings.append(dict(route=route,image_id=image_id,scope=domain,wall_s=elapsed))
                for i,box in enumerate(f['boxes']):
                    y0,x0,y1,x1=box
                    if domain=='ROI':
                        tick=time.monotonic();objs,_,labs=engine.detect(f['image'][y0:y1,x0:x1])
                        elapsed=time.monotonic()-tick
                        timings.append(dict(route=route,image_id=image_id,roi_index=i,scope=domain,wall_s=elapsed))
                    else:
                        labs=labels[y0:y1,x0:x1]
                        objs=[dict(o,xy_local=[o['xy_local'][0]-x0,o['xy_local'][1]-y0]) for o in objects
                              if x0<=o['xy_local'][0]<x1 and y0<=o['xy_local'][1]<y1]
                    metrics=evaluate_instances(full[y0:y1,x0:x1],labs,objs,
                        (f['masks']['UNCERTAIN'][y0:y1,x0:x1]>0)|(f['masks']['OUT_OF_FOCUS'][y0:y1,x0:x1]>0),counts)
                    key=f'{route}_{domain}_{image_id}_roi{i:02d}'
                    result=dict(model=route,split='test',domain=domain,image_id=image_id,roi_index=i,box_yxyx=box,
                                metrics=metrics,objects=objs,detect_wall_s=elapsed,artifact=key,
                                annotation_revision=f['row']['annotation']['revision'])
                    np.savez_compressed(raw/(key+'.npz'),labels=labs,truth=full[y0:y1,x0:x1])
                    write(raw/(key+'.json'),result);rows.append(result)
                write(out/'rows.partial.json',rows)
        print('Completed',route,flush=True)
    metrics={f'{r}/{d}':aggregate([a for a in rows if a['model']==r and a['domain']==d])
             for r in configs for d in ['ROI','FULL_FRAME']}
    per_frame={f'{r}/{image_id}/{d}':aggregate([a for a in rows if (a['model'],a['image_id'],a['domain'])==(r,image_id,d)])
               for r in configs for image_id in selected_ids for d in ['ROI','FULL_FRAME']}
    for route in configs:
        assert sha(BASE/f'selected_{route}.json')==prereg['config_sha256'][route]
    assert sha(ROOT/'models/cellpose/cyto3')==EXPECTED_CYTO3
    write(out/'all_rows.json',rows);write(out/'timings.json',timings)
    result=dict(status='COMPLETED',metrics=metrics,per_frame=per_frame,coverage=coverage,
                source_frame_count=2,roi_count=sum(len(f['boxes']) for f in frames),
                reviewed_single_masks=sum(c['reviewed_counts']['SINGLE_BUBBLE'] for c in coverage),
                unreviewed_saved_single_masks=sum(c['unreviewed_saved_single_masks'] for c in coverage),
                wall_s=time.monotonic()-started,configs_unchanged=True,model_sha256_unchanged=True,
                no_test_training=True,no_test_tuning=True,no_test_model_selection=True,
                BUBBLE_DETECTION_VALIDATED=False,SHORT_TRACKING_VALIDATED=False,GROUP_SEGMENTATION_VALIDATED=False,
                group_shape=None,control_permission=False,training_ready=False)
    write(out/'summary.json',result)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--freeze',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    def timeout(*_):raise TimeoutError('HELDOUT_WALL_BUDGET_180S')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
    try:evaluate(args.freeze.resolve(),args.out.resolve())
    except Exception:
        if args.out.exists():write(args.out/'FAILED.json',{'error':traceback.format_exc()})
        raise
    finally:signal.alarm(0)


if __name__=='__main__':main()

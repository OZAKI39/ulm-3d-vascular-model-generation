"""Original SAM: cancellable process, reusable weights and bounded feature LRU.

JSON-lines server plus a one-request replay CLI. Original masks retained;
point-instance mode explicitly extracts the clicked pixel's connected component.
"""
from pathlib import Path
import sys, json, time, traceback, select
from collections import OrderedDict
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))


def emit(**message): print(json.dumps(message, ensure_ascii=False), flush=True)


class SamWorker:
    def __init__(self): self.predictor = None; self.cache = OrderedDict()

    def run(self, req):
        if req.get('mode')=='point_instances':return self.run_instances(req)
        start = time.monotonic(); token = req['token']
        def status(text): emit(token=token, status=text)
        import numpy as np, imageio.v3 as iio, torch
        from micro_sam import util
        from sam_prompt_geometry import prompt_region, proposal_checks, single_bubble_checks
        torch.set_num_threads(2); cold = self.predictor is None
        if cold:
            status('首次加载原版 SAM ViT-B；随后复用模型，可取消')
            self.predictor = util.get_sam_model(model_type='vit_b', device='cuda' if torch.cuda.is_available() else 'cpu',
                                               checkpoint_path=str(ROOT / 'models/sam_vit_b_01ec64.pth'))
        predictor = self.predictor
        image_path = Path(req['image']); image = iio.imread(image_path)
        positive = np.asarray(req['positive_yx'], float).reshape(-1, 2)
        negative = np.asarray(req['negative_yx'], float).reshape(-1, 2)
        local = req.get('mode') == 'local_box'
        if local and req.get('object_class')=='SINGLE_BUBBLE':
            width=req.get('target_width_px',16)
            if len(positive)>1 and np.max(np.ptp(positive,axis=0))>width:
                raise ValueError('这些绿点相距较远，请使用“每点一个单泡”模式，不能合成一个单泡')
            b=req.get('box_yxyx')
            if b is not None and max(b[2]-b[0],b[3]-b[1])>2*width:
                raise ValueError('单个对象的框太大；多泡请切换“每点一个单泡”，或重新紧框一个泡')
        if local:
            crop, box = prompt_region(image.shape, positive, negative, req.get('box_yxyx'),
                                      req.get('target_width_px', 32), req.get('context_width_px', 128))
        else:
            crop = np.array([0, 0, *image.shape[:2]]); box = None
        y0, x0, y1, x1 = crop; offset = np.array([y0, x0])
        key = (str(image_path.resolve()), image_path.stat().st_mtime_ns, tuple(crop))
        cache_hit = key in self.cache
        status('复用局部图像特征' if cache_hit else '计算局部图像特征' if local else '计算整图特征（旧版对照）')
        if cache_hit: embeddings = self.cache.pop(key)
        else: embeddings = util.precompute_image_embeddings(predictor, image[y0:y1, x0:x1], ndim=2, verbose=False)
        self.cache[key] = embeddings
        while len(self.cache) > 16: self.cache.popitem(last=False)
        util.set_precomputed(predictor, embeddings)
        pts = np.concatenate([positive, negative], axis=0) - offset
        labels = np.r_[np.ones(len(positive), np.int32), np.zeros(len(negative), np.int32)]
        local_box = box - np.r_[offset, offset] if box is not None else None
        status('生成多个候选；只检查提示约束，不自动确认标签')
        masks, scores, _ = predictor.predict(point_coords=pts[:, ::-1] if len(pts) else None,
                                            point_labels=labels if len(pts) else None,
                                            box=local_box[[1, 0, 3, 2]] if local_box is not None else None,
                                            multimask_output=local)
        info = [proposal_checks(m, positive-offset, negative-offset, local_box) if local else
                {'allowed': True, 'reasons': [], 'area_px': int(m.sum())} for m in masks]
        if local and req.get('object_class')=='SINGLE_BUBBLE':
            info=[single_bubble_checks(m,positive-offset,negative-offset,local_box,req.get('target_width_px',16)) for m in masks]
        order = sorted(range(len(masks)), key=lambda i: (not info[i]['allowed'], -float(scores[i])))
        full = np.zeros((len(masks), *image.shape[:2]), bool)
        for n, i in enumerate(order): full[n, y0:y1, x0:x1] = masks[i]
        info = [{**info[i], 'model_score': float(scores[i])} for i in order]
        result = {'wall_s': time.monotonic()-start, 'cold_model': cold, 'embedding_cache_hit': cache_hit,
                  'crop_yxyx': crop.tolist(), 'box_yxyx': box.tolist() if box is not None else None,
                  'proposals': info, 'mask_geometry_modified': False, 'mode': req.get('mode', 'legacy_full_frame')}
        np.savez_compressed(req['output'], mask=full[0], masks=full, score=scores[order],
                            allowed=np.array([v['allowed'] for v in info]), details=json.dumps(result, ensure_ascii=False),
                            request_token=token, image_id=req['image_id'])
        emit(token=token, event='result', output=req['output'], **result)

    def run_instances(self,req):
        import numpy as np,imageio.v3 as iio
        from sam_prompt_geometry import point_instance_requests,point_component,single_bubble_checks
        start=time.monotonic();image=iio.imread(req['image']);cold=self.predictor is None
        plans=point_instance_requests(image.shape,req['positive_yx'],req['negative_yx'],req.get('box_yxyx'),req.get('target_width_px',16))
        folder=Path(req['output']).parent/'per_point';folder.mkdir(exist_ok=True)
        arrays={};targets=[]
        for plan in plans:
            i=plan['index'];emit(token=req['token'],status=f'逐泡分割 {i+1}/{len(plans)}：每个绿点独立计算')
            child={**req,'mode':'local_box','object_class':'SINGLE_BUBBLE','box_yxyx':plan['box_yxyx'],
                   'positive_yx':plan['positive_yx'],'negative_yx':plan['negative_yx'],
                   'token':f'{req["token"]}_point_{i}','output':str(folder/f'{i:03d}.npz')}
            (folder/f'{i:03d}_request.json').write_text(json.dumps(child,indent=2))
            self.run(child)
            with np.load(child['output'],allow_pickle=False) as z:
                details=json.loads(str(z['details']));y0,x0,y1,x1=details['crop_yxyx']
                raw=z['masks'][:,y0:y1,x0:x1].copy()
            local_positive=np.asarray(plan['positive_yx'])-[y0,x0]
            local_negative=np.asarray(plan['negative_yx']).reshape(-1,2)-[y0,x0]
            local_box=np.asarray(plan['box_yxyx'])-[y0,x0,y0,x0]
            processed=[];checks=[]
            for mask,old_info in zip(raw,details['proposals']):
                instance,extraction=point_component(mask,local_positive[0]);processed.append(instance)
                info=single_bubble_checks(instance,local_positive,local_negative,local_box,req.get('target_width_px',16))
                checks.append({**info,**extraction,'model_score':old_info['model_score']})
            order=sorted(range(len(checks)),key=lambda j:(not checks[j]['allowed'],-checks[j]['model_score']))
            arrays[f'target_{i:03d}']=np.asarray(processed)[order]
            arrays[f'raw_target_{i:03d}']=raw[order]
            details['proposals']=[checks[j] for j in order]
            details['mask_geometry_modified']=any(p['removed_other_component_area_px']>0 for p in checks)
            details['postprocessing']='clicked_pixel_8_connected_component; raw masks retained'
            if plan['duplicate_point_indices']:
                for p in details['proposals']:
                    p['allowed']=False;p['reasons'].append('附近有重复绿点；同一个泡只点一次')
            targets.append({**plan,**details,'array_key':f'target_{i:03d}'})
        result={'mode':'point_instances','wall_s':time.monotonic()-start,'cold_model':cold,
                'target_count':len(targets),'admissible_target_count':sum(any(p['allowed'] for p in t['proposals']) for t in targets),
                'targets':targets,'image_shape_yx':list(image.shape[:2]),'roi_box_yxyx':req.get('box_yxyx'),
                'box_semantics':'selection_ROI_only; independent target boxes per point',
                'mask_geometry_modified':any(p['removed_other_component_area_px']>0 for t in targets for p in t['proposals']),
                'postprocessing':'clicked_pixel_component_only; raw masks retained; no hole fill/smoothing/box clipping',
                'embedding_cache_hit':all(t['embedding_cache_hit'] for t in targets)}
        np.savez_compressed(req['output'],**arrays,details=json.dumps(result,ensure_ascii=False),
                            request_token=req['token'],image_id=req['image_id'])
        emit(token=req['token'],event='result',output=req['output'],**result)


def main():
    worker = SamWorker()
    if sys.argv[1:] == ['--server']:
        emit(event='ready')
        while select.select([sys.stdin], [], [], 120)[0]:
            line = sys.stdin.readline()
            if not line: break
            req = None
            try:
                req = json.loads(line); worker.run(req)
            except Exception:
                emit(event='error', token=req.get('token') if req else None, error=traceback.format_exc())
        emit(event='idle_exit')
    else: worker.run(json.loads(Path(sys.argv[1]).read_text()))

if __name__ == '__main__': main()

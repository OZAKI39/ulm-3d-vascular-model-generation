"""Original SAM: cancellable process, reusable weights and four-crop feature LRU.

JSON-lines server plus a one-request replay CLI. No SAM2 or altered silhouettes.
"""
from pathlib import Path
import sys, json, time, traceback, select
from collections import OrderedDict
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))


def emit(**message): print(json.dumps(message, ensure_ascii=False), flush=True)


class SamWorker:
    def __init__(self): self.predictor = None; self.cache = OrderedDict()

    def run(self, req):
        start = time.monotonic(); token = req['token']
        def status(text): emit(token=token, status=text)
        import numpy as np, imageio.v3 as iio, torch
        from micro_sam import util
        from sam_prompt_geometry import prompt_region, proposal_checks
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
        while len(self.cache) > 4: self.cache.popitem(last=False)
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

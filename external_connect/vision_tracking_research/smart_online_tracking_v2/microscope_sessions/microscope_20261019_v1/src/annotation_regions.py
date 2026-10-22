"""Explicit spatial review coverage. Unreviewed pixels are never background truth."""
import hashlib,json
import numpy as np


def normalize_box(box,shape):
    a=np.asarray(box,dtype=float)
    if a.shape!=(4,) or not np.isfinite(a).all() or not np.equal(a,np.round(a)).all():
        raise ValueError('区域必须是四个整数像素坐标 y0,x0,y1,x1')
    y0,x0,y1,x1=map(int,a)
    if not (0<=y0<y1<=shape[0] and 0<=x0<x1<=shape[1]):
        raise ValueError('审核区域必须完整位于图像内')
    return [y0,x0,y1,x1]


def region_digest(masks,box):
    shape=next(iter(masks.values())).shape
    y0,x0,y1,x1=normalize_box(box,shape)
    h=hashlib.sha256(json.dumps([y0,x0,y1,x1]).encode())
    for name in sorted(masks):
        h.update(name.encode());h.update(np.asarray(masks[name][y0:y1,x0:x1],dtype='<u4').tobytes())
    return h.hexdigest()


def valid_regions(meta,masks):
    """Fail closed on changed or malformed region proofs; preserve other regions."""
    result=[]
    for r in (meta or {}).get('review_regions',[]):
        try:
            if not r.get('reviewer','').strip() or not r.get('confirmed',False):continue
            box=normalize_box(r['box_yxyx'],next(iter(masks.values())).shape)
            if region_digest(masks,box)==r['mask_sha256']:result.append(dict(r,box_yxyx=box))
        except (KeyError,ValueError,TypeError):continue
    return result


def reviewed_boxes(meta,masks):
    shape=next(iter(masks.values())).shape
    if not meta:return []
    if meta.get('expert_confirmed') and meta.get('status')=='EXPERT_CONFIRMED' and meta.get('review_scope','FULL_FRAME')=='FULL_FRAME':
        return [[0,0,*shape]]
    return [r['box_yxyx'] for r in valid_regions(meta,masks)]


def review_identity(meta):
    return {'scope':meta.get('review_scope','FULL_FRAME'),'full_confirmed':meta.get('expert_confirmed',False),
            'regions':[{k:r.get(k) for k in ('id','box_yxyx','mask_sha256','reviewer','confirmed')}
                       for r in meta.get('review_regions',[])]}


def reviewed_counts(meta,masks):
    coverage=np.zeros(next(iter(masks.values())).shape,bool)
    for y0,x0,y1,x1 in reviewed_boxes(meta,masks):coverage[y0:y1,x0:x1]=True
    return {k:int(np.count_nonzero(np.unique(a[coverage]))) for k,a in masks.items()}

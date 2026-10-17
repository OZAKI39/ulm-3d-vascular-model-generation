"""Auditable instance edits. No automatic conversion of predictions to reviewed GT."""
import numpy as np
from scipy import ndimage as ndi
from skimage.segmentation import watershed


def merge_proposals(masks, proposals, box):
    """Add whole, non-overlapping instances, preserving every existing pixel."""
    y0,x0,y1,x1=box
    if proposals.shape!=(y1-y0,x1-x0):raise ValueError('Candidate shape mismatch')
    target=masks['SINGLE_BUBBLE'].copy()
    occupied=np.logical_or.reduce([v[y0:y1,x0:x1]>0 for v in masks.values()])
    local=target[y0:y1,x0:x1];next_id=int(target.max())+1;accepted=[];skipped=[]
    for ident in np.unique(proposals):
        if not ident:continue
        mask=proposals==ident
        if np.any(mask&occupied):skipped.append(int(ident));continue
        local[mask]=next_id;occupied|=mask
        accepted.append({'candidate_id':int(ident),'instance_id':next_id});next_id+=1
    return target,accepted,skipped


def split_instance(labels, ident, points):
    """Two human seeds split the selected mask; its exact pixel union is preserved."""
    mask=labels==ident
    if not ident or not mask.any():raise ValueError('请先选中一个已有对象')
    if ndi.label(mask)[1]!=1:raise ValueError('此编号含多个不相连区域，请先用画笔修正或分别编号')
    points=np.asarray(points,float)
    if points.shape!=(2,2) or not np.isfinite(points).all():raise ValueError('请在两个泡内部各放一个点，共两个点')
    pts=np.floor(points+.5).astype(int)
    if np.array_equal(pts[0],pts[1]):raise ValueError('两个点不能重合')
    markers=np.zeros(labels.shape,np.int32)
    for i,(y,x) in enumerate(pts,1):
        if not (0<=y<labels.shape[0] and 0<=x<labels.shape[1] and mask[y,x]):
            raise ValueError('两个分离点都必须在选中对象内部')
        markers[y,x]=i
    parts=watershed(-ndi.distance_transform_edt(mask),markers,mask=mask)
    if not np.array_equal(parts>0,mask) or not all(np.any(parts==i) for i in (1,2)):
        raise ValueError('无法完整分开，请重新放置两个点')
    result=labels.copy();new_id=int(labels.max())+1
    result[parts==2]=new_id
    return result,new_id

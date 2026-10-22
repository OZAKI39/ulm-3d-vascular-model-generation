"""Original-pixel prompts, explicit seed-component extraction and proposal guards."""
import numpy as np


def point_component(mask,point_yx):
    """Select the raw 8-connected component containing the clicked pixel.

    Does not fill holes, smooth, clip to a box, or relocate a missed prompt.
    Caller retains the complete raw SAM mask alongside this instance mask.
    """
    from scipy.ndimage import label
    raw=np.asarray(mask,bool);ij=np.clip(np.floor(np.asarray(point_yx)+.5).astype(int),[0,0],np.array(raw.shape)-1)
    components,count=label(raw,structure=np.ones((3,3),np.uint8));component=int(components[tuple(ij)])
    instance=components==component if component else np.zeros_like(raw)
    return instance,{'raw_area_px':int(raw.sum()),'raw_component_count':int(count),
                     'removed_other_component_area_px':int(raw.sum()-instance.sum()),
                     'postprocessing':'clicked_pixel_8_connected_component','seed_in_raw_mask':bool(component)}


def prompt_region(shape, positive_yx, negative_yx, box_yxyx=None, target_width=32, context_width=128):
    h, w = shape[:2]
    positive = np.asarray(positive_yx, dtype=float).reshape(-1, 2)
    negative = np.asarray(negative_yx, dtype=float).reshape(-1, 2)
    for points in (positive, negative):
        if not np.isfinite(points).all() or np.any(points < 0) or np.any(points >= [h, w]):
            raise ValueError('提示点必须位于原图内部')
    if box_yxyx is None:
        if not len(positive): raise ValueError('请框选目标，或在一个对象内部放正点')
        lo = positive.min(axis=0) - target_width / 2
        hi = positive.max(axis=0) + target_width / 2
    else:
        b = np.asarray(box_yxyx, dtype=float)
        if b.shape != (4,) or not np.isfinite(b).all(): raise ValueError('目标框无效')
        lo, hi = b[:2], b[2:]
    lo = np.maximum(0, lo); hi = np.minimum([h, w], hi)
    if np.any(hi - lo < 3): raise ValueError('目标框太小；请框住一个完整对象')
    if len(positive) and (np.any(positive < lo) or np.any(positive >= hi)):
        raise ValueError('正点在目标框外；请调整框或清空旧点')
    center = (lo + hi) / 2; side = np.maximum(context_width, 2 * (hi - lo))
    crop_lo = np.maximum(0, np.floor(center - side / 2)).astype(int)
    crop_hi = np.minimum([h, w], np.ceil(center + side / 2)).astype(int)
    if np.any(crop_hi - crop_lo > 1024): raise ValueError('目标范围过大；请只框选当前一个对象')
    if len(negative) and (np.any(negative < crop_lo) or np.any(negative >= crop_hi)):
        raise ValueError('负点离目标太远；清空旧点后在紧邻误选区域放点')
    return np.r_[crop_lo, crop_hi], np.r_[lo, hi]


def proposal_checks(mask, positive_yx, negative_yx, box_yxyx):
    mask = np.asarray(mask, bool); reasons = []; area = int(mask.sum())
    positive = np.asarray(positive_yx).reshape(-1, 2); negative = np.asarray(negative_yx).reshape(-1, 2)
    if not area: reasons.append('空候选')
    for points, desired, text in [(positive, True, '漏掉正点'), (negative, False, '包含负点')]:
        if len(points):
            ij = np.clip(np.floor(points + .5).astype(int), [0, 0], np.array(mask.shape) - 1)
            if np.any(mask[ij[:, 0], ij[:, 1]] != desired): reasons.append(text)
    if area > .70 * mask.size: reasons.append('覆盖局部视野过多，疑似通道')
    if mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any():
        reasons.append('触及局部裁剪边缘，请重新框选')
    y0, x0, y1, x1 = np.asarray(box_yxyx, float)
    y0, x0 = np.floor([y0, x0]).astype(int); y1, x1 = np.ceil([y1, x1]).astype(int)
    inside = int(mask[max(y0, 0):min(y1, mask.shape[0]), max(x0, 0):min(x1, mask.shape[1])].sum())
    if area and inside / area < .8: reasons.append('大部分轮廓超出目标框')
    return {'allowed': not reasons, 'reasons': reasons, 'area_px': area,
            'crop_fraction': area / mask.size, 'inside_box_fraction': inside / max(area, 1)}


def single_bubble_checks(mask, positive_yx, negative_yx, box_yxyx, target_width):
    """A user-adjustable pixel-size guard, not a biological classifier."""
    result=proposal_checks(mask,positive_yx,negative_yx,box_yxyx)
    yy,xx=np.where(mask)
    span=[int(np.ptp(yy)+1),int(np.ptp(xx)+1)] if len(yy) else [0,0]
    if result['area_px']>1.5*target_width**2 or max(span)>1.5*target_width:
        result['reasons'].append('超过当前单泡局部尺寸；不会按大框放宽限制')
    result.update(allowed=not result['reasons'],span_yx_px=span,target_width_px=target_width)
    return result


def point_instance_requests(shape, positive_yx, negative_yx, roi_box, target_width, context_width=128):
    """Each positive point is a separate instance; a large box is only an ROI.

    Nearby other positive points become exclusion prompts for this instance.
    Negative points outside a local crop are recorded as not applicable.
    """
    positive=np.asarray(positive_yx,float).reshape(-1,2)
    negative=np.asarray(negative_yx,float).reshape(-1,2)
    if not 1<=len(positive)<=64:raise ValueError('单泡点选需要 1–64 个绿点，每点对应一个泡')
    for points in [positive,negative]:
        if not np.isfinite(points).all() or np.any(points<0) or np.any(points>=np.array(shape[:2])):
            raise ValueError('提示点必须在原图内')
    if roi_box is not None:
        box=np.asarray(roi_box,float)
        if box.shape!=(4,) or not np.isfinite(box).all() or np.any(box[2:]<=box[:2]):raise ValueError('选区框无效')
        if np.any(positive<box[:2]) or np.any(positive>=box[2:]):raise ValueError('有绿点在选区框外；请扩大或清除选区框')
    plans=[]
    for i,point in enumerate(positive):
        crop,box=prompt_region(shape,[point],[],target_width=target_width,context_width=context_width)
        nearby=np.all((negative>=crop[:2])&(negative<crop[2:]),axis=1)
        others=np.arange(len(positive))!=i
        adjacent=others & np.all((positive>=crop[:2])&(positive<crop[2:]),axis=1)
        duplicates=np.flatnonzero(others & (np.linalg.norm(positive-point,axis=1)<2))
        plans.append({'index':i,'positive_yx':[point.tolist()],
                      'negative_yx':np.concatenate([negative[nearby],positive[adjacent]],axis=0).tolist(),
                      'user_negative_indices':np.flatnonzero(nearby).tolist(),
                      'excluded_other_point_indices':np.flatnonzero(adjacent).tolist(),
                      'duplicate_point_indices':duplicates.tolist(),
                      'crop_yxyx':crop.tolist(),'box_yxyx':box.tolist()})
    return plans

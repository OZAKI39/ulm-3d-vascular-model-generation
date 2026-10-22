"""Original-pixel prompt geometry and proposal guards; never reshape SAM masks."""
import numpy as np


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

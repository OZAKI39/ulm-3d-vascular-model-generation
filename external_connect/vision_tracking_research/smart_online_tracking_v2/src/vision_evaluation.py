"""One-to-one instance evaluation only inside explicitly reviewed original pixels."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def assignment(score, eligible):
    if not score.size:return []
    # Cardinality takes priority; invalid pairs cannot displace a valid match.
    reward=eligible*(max(score.shape)+1)+np.where(eligible,score,0)
    rr,cc=linear_sum_assignment(reward,maximize=True)
    return [(int(r),int(c)) for r,c in zip(rr,cc) if eligible[r,c]]


def rates(tp,fp,fn):
    return {'tp':tp,'fp':fp,'fn':fn,'precision':tp/(tp+fp) if tp+fp else None,
            'recall':tp/(tp+fn) if tp+fn else None,'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None}


def evaluate_instances(gt,pred,objects,ignore=None,full_gt_counts=None,margin=6,radius=6,iou_threshold=.5):
    """Objects use local xy centers. Partial/border/ignore predictions are counted separately.

    All accepted prediction objects have area <=450 under the inherited candidate
    classifier. Large regions remain unverified group candidates, not single TPs.
    """
    h,w=gt.shape;ignore=np.zeros(gt.shape,bool) if ignore is None else np.asarray(ignore,bool)
    gids=[];gcenters=[];excluded_gt=[]
    for label in np.unique(gt):
        if not label:continue
        y,x=np.where(gt==label);center=np.array([x.mean(),y.mean()])
        complete=(full_gt_counts is None or len(x)==full_gt_counts[int(label)])
        if not complete or not(margin<=center[0]<w-margin and margin<=center[1]<h-margin):
            excluded_gt.append(int(label));continue
        gids.append(int(label));gcenters.append(center)
    boundary_truth=np.isin(gt,excluded_gt);ignored=[];pids=[];pcenters=[];large=[]
    for obj in objects:
        label=int(obj['label']);mask=pred==label
        if not np.any(mask):continue
        if obj['area']>450:large.append(label);continue
        x,y=obj['xy_local']
        reason=None
        if not(margin<=x<w-margin and margin<=y<h-margin):reason='BORDER_CENTER'
        elif np.count_nonzero(mask & (ignore|boundary_truth))/np.count_nonzero(mask)>=.5:reason='IGNORE_OR_TRUNCATED_TRUTH'
        if reason:ignored.append({'label':label,'reason':reason});continue
        pids.append(label);pcenters.append([x,y])
    g=np.asarray(gcenters).reshape(-1,2);p=np.asarray(pcenters).reshape(-1,2)
    distance=np.linalg.norm(g[:,None,:]-p[None,:,:],axis=2)
    center_pairs=assignment(1-np.minimum(distance/radius,1),distance<=radius)
    iou=np.zeros((len(gids),len(pids)))
    for i,gid in enumerate(gids):
        a=gt==gid
        for j,pid in enumerate(pids):
            b=pred==pid;iou[i,j]=np.count_nonzero(a&b)/np.count_nonzero(a|b)
    mask_pairs=assignment(iou,iou>=iou_threshold)
    matched_g={i for i,_ in center_pairs};matched_p={j for _,j in center_pairs}
    return {'center':rates(len(center_pairs),len(pids)-len(center_pairs),len(gids)-len(center_pairs)),
            'mask':rates(len(mask_pairs),len(pids)-len(mask_pairs),len(gids)-len(mask_pairs)),
            'center_errors_px':[float(distance[i,j]) for i,j in center_pairs],
            'matched_ious':[float(iou[i,j]) for i,j in mask_pairs],
            'matches':[{'gt':gids[i],'pred':pids[j],'distance_px':float(distance[i,j])} for i,j in center_pairs],
            'false_positive_ids':[pids[i] for i in range(len(pids)) if i not in matched_p],
            'false_negative_ids':[gids[i] for i in range(len(gids)) if i not in matched_g],
            'eligible_gt':len(gids),'eligible_predictions':len(pids),'excluded_gt_ids':excluded_gt,
            'ignored_predictions':ignored,'large_unverified_regions':large,
            'ignore_pixels':int(ignore.sum()),'reviewed_pixels':int(gt.size)}


def aggregate(rows):
    result={'roi_count':len(rows),'source_frame_count':len({r['image_id'] for r in rows})}
    for kind in ('center','mask'):
        counts=[sum(r['metrics'][kind][k] for r in rows) for k in ('tp','fp','fn')]
        result[kind]=rates(*counts)
    errors=[v for r in rows for v in r['metrics']['center_errors_px']]
    result['center_error_px_mean']=float(np.mean(errors)) if errors else None
    result['center_error_px_p95']=float(np.percentile(errors,95)) if errors else None
    result['fp_per_reviewed_roi']=result['center']['fp']/len(rows) if rows else None
    result['fn_per_reviewed_roi']=result['center']['fn']/len(rows) if rows else None
    result['excluded_gt']=sum(len(r['metrics']['excluded_gt_ids']) for r in rows)
    result['ignored_predictions']=sum(len(r['metrics']['ignored_predictions']) for r in rows)
    result['large_unverified_regions']=sum(len(r['metrics']['large_unverified_regions']) for r in rows)
    result['ignore_pixels']=sum(r['metrics']['ignore_pixels'] for r in rows)
    result['uncertain_rejection_rate']=None if result['ignore_pixels']==0 else 'See per-ROI exclusions'
    return result

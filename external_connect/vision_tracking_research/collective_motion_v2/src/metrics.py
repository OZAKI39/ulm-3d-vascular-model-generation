"""Explicit-denominator scientific metrics; unchanged F10U tolerances."""
import numpy as np
from vcommon import LIMITS,NAMES

def rate(num,den):return float(num/den) if den else None

def summarize(mean,truth,label,accepted,predicted_class,score,eligible,metadata,sigma=None,quantiles=None):
    mean=np.asarray(mean,float);truth=np.asarray(truth,float);label=np.asarray(label);accepted=np.asarray(accepted,bool)
    finite=np.isfinite(mean).all(1);global_truth=label==0;err=mean-truth;abs_err=np.abs(err)
    bad=~global_truth|~finite|np.any(abs_err>LIMITS,axis=1);valid_truth=global_truth&finite
    out=dict(n_total=len(label),n_global_truth=int(global_truth.sum()),n_non_affine=int((label==1).sum()),n_not_observable=int((label==2).sum()),n_finite_predictions=int(finite.sum()),
        n_accepted=int(accepted.sum()),acceptance_rate=rate(accepted.sum(),len(label)),n_accepted_failed=int((accepted&bad).sum()),accepted_failure_fraction=rate((accepted&bad).sum(),accepted.sum()),
        n_accepted_global_truth=int((accepted&global_truth).sum()),n_accepted_truth_over_tolerance=int((accepted&global_truth&bad).sum()),accepted_truth_failure_fraction=rate((accepted&global_truth&bad).sum(),(accepted&global_truth).sum()),
        n_accurate_accepted=int((accepted&~bad).sum()),accurate_accepted_fraction=rate((accepted&~bad).sum(),len(label)),n_correct_rejection_non_global=int((~accepted&~global_truth).sum()),correct_rejection_non_global_fraction=rate((~accepted&~global_truth).sum(),(~global_truth).sum()),
        n_rejected_global=int((~accepted&global_truth).sum()),over_rejection_fraction_global=rate((~accepted&global_truth).sum(),global_truth.sum()),n_global_with_any_error=int((bad&global_truth).sum()),all_global_failure_fraction=rate((bad&global_truth).sum(),global_truth.sum()),
        n_mae_denominator=int(valid_truth.sum()),mae_all_finite_global=np.mean(abs_err[valid_truth],axis=0).tolist() if valid_truth.any() else [None]*6,
        bias_all_finite_global=np.mean(err[valid_truth],axis=0).tolist() if valid_truth.any() else [None]*6,
        mae_accepted_global=np.mean(abs_err[valid_truth&accepted],axis=0).tolist() if (valid_truth&accepted).any() else [None]*6,
        per_parameter_over_limit_counts=[int(np.sum(global_truth&(~finite|(abs_err[:,j]>LIMITS[j])))) for j in range(6)])
    signs=[]
    for j in range(6):
        mask=valid_truth&(np.abs(truth[:,j])>LIMITS[j]*.25);right=(np.sign(truth[:,j])==np.sign(mean[:,j]))&mask
        signs.append(dict(n=int(mask.sum()),correct=int(right.sum()),accuracy=rate(right.sum(),mask.sum()),minimum_truth_magnitude=float(LIMITS[j]*.25)))
    out['sign_accuracy']=signs
    confusion=np.zeros((3,3),int)
    for a,b in zip(label,predicted_class):confusion[a,b]+=1
    out['class_confusion_true_rows_predicted_columns']=confusion.tolist()
    order=np.argsort(score);order=order[np.asarray(eligible,bool)[order]&finite[order]&np.isfinite(score[order])]
    risk=np.cumsum(bad[order])/np.arange(1,len(order)+1)
    out['risk_coverage']=dict(n_total=len(label),n_rankable=len(order),coverage=(np.arange(1,len(order)+1)/len(label)).tolist(),risk=risk.tolist(),ordered_indices=order.tolist(),scope='vary score threshold among algebra-valid predicted-global candidates; operational validation threshold reported separately')
    out['by_case']={};out['by_source']={}
    for grouping,key in [('by_case','case'),('by_source','source_frame')]:
        for value in sorted(set(m[key] for m in metadata)):
            m=np.array([z[key]==value for z in metadata]);g=m&valid_truth
            out[grouping][str(value)]=dict(n=int(m.sum()),n_truth=int((m&global_truth).sum()),n_accepted=int((m&accepted).sum()),n_failed_accepted=int((m&accepted&bad).sum()),accepted_failure_fraction=rate((m&accepted&bad).sum(),(m&accepted).sum()),mae=np.mean(abs_err[g],axis=0).tolist() if g.any() else [None]*6)
    if sigma is not None:
        out['intervals']={}
        for level,q in (quantiles or {}).items():
            q=np.asarray(q,float);width=np.asarray(sigma)*q
            coverage=abs_err<=width
            out['intervals'][level]=dict(n=int(valid_truth.sum()),covered_counts=coverage[valid_truth].sum(0).tolist(),coverage=coverage[valid_truth].mean(0).tolist(),joint_coverage=float(coverage[valid_truth].all(1).mean()),mean_width=(2*width[valid_truth]).mean(0).tolist(),nominal_level=float(level),real_domain_calibrated=False,by_case={c:dict(n=int((valid_truth&np.array([m['case']==c for m in metadata])).sum()),coverage=coverage[valid_truth&np.array([m['case']==c for m in metadata])].mean(0).tolist()) for c in sorted(set(m['case'] for m in metadata if m['class_id']==0))})
    sources=np.array([m['source_frame'] for m in metadata]);unique=np.unique(sources);rng=np.random.default_rng(29171);boot=[]
    for _ in range(500):
        ids=np.concatenate([np.flatnonzero(sources==s) for s in rng.choice(unique,len(unique),replace=True)])
        boot.append([np.mean(accepted[ids]),rate(np.sum(accepted[ids]&bad[ids]),np.sum(accepted[ids]))])
    out['source_block_bootstrap']=dict(n_source_frames=len(unique),independent_recordings=1,claim='descriptive resampling of four source frames; NOT independent-experiment confidence',acceptance_interval95=np.percentile(np.array(boot,dtype=float)[:,0],[2.5,97.5]).tolist(),accepted_risk_interval95=np.nanpercentile(np.array(boot,dtype=float)[:,1],[2.5,97.5]).tolist() if any(b[1] is not None for b in boot) else None)
    return out

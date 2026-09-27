"""Held-out semantic and binary ROI scores; UNKNOWN remains a false negative."""
import numpy as np
from scipy.spatial import cKDTree

def scores(truth,prediction):
    truth,prediction=np.asarray(truth),np.asarray(prediction)
    f1=[]
    for label in [1,2,3]:
        tp=np.sum((truth==label)&(prediction==label));fp=np.sum((truth!=label)&(prediction==label));fn=np.sum((truth==label)&(prediction!=label))
        f1.append(float(2*tp/(2*tp+fp+fn)) if 2*tp+fp+fn else 0.)
    return dict(M1_f1=f1[0],M2_f1=f1[1],M3_f1=f1[2],macro_f1=float(np.mean(f1)),UNKNOWN_fraction=float(np.mean(prediction==0)))


def boundary_points(points,labels,pair):
    """Evaluation only: label interfaces on a local geometric k-neighbour graph.

    No class labels or boundary coordinates enter registration or transport.
    The same local geometric graph is used for true and predicted interfaces.
    """
    count=min(7,len(points));d,index=cKDTree(points).query(points,k=count)
    local_scale=np.median(d[:,1])
    edges=set()
    for i in range(len(points)):
        for dist,j in zip(d[i,1:],index[i,1:]):
            if dist<=3*local_scale and {int(labels[i]),int(labels[j])}==set(pair):edges.add(tuple(sorted([i,int(j)])))
    return np.array([(points[i]+points[j])/2 for i,j in sorted(edges)]).reshape(-1,3)


def boundary_error(points,truth,prediction,pair):
    actual=boundary_points(points,truth,pair);predicted=boundary_points(points,prediction,pair)
    if not len(actual) or not len(predicted):return None
    return float(np.median(np.r_[cKDTree(actual).query(predicted)[0],cKDTree(predicted).query(actual)[0]]))


def collapse_labels_to_mevo(labels, *, truth=False):
    """0=UNKNOWN, 1=M1, 2=MeVO; preserve the existing rejection decision."""
    labels=np.asarray(labels)
    allowed=[1,2,3] if truth else [0,1,2,3]
    if labels.ndim!=1 or not np.isin(labels,allowed).all():
        raise ValueError('Expected a vector of native semantic labels; truth cannot be UNKNOWN')
    return np.where(labels==3,2,labels).astype(np.int8)


def evaluate_mevo_binary(truth, prediction):
    actual=collapse_labels_to_mevo(truth,truth=True)
    predicted=collapse_labels_to_mevo(prediction)
    if actual.shape!=predicted.shape or not len(actual):
        raise ValueError('Truth and prediction must cover the same nonempty sample set')
    # UNKNOWN contributes to the true class FN through the full row total.
    matrix=np.array([[np.sum((actual==a)&(predicted==b)) for b in [1,2,0]] for a in [1,2]])
    divide=lambda a,b:float(a/b) if b else 0.
    result={}
    for i,name in enumerate(['M1','MeVO']):
        tp=matrix[i,i];fn=matrix[i].sum()-tp;fp=matrix[:,i].sum()-tp
        result.update({name+'_precision':divide(tp,tp+fp),name+'_recall':divide(tp,tp+fn),
                       name+'_f1':divide(2*tp,2*tp+fp+fn)})
    result.update(binary_macro_f1=(result['M1_f1']+result['MeVO_f1'])/2,
        false_inclusion=divide(matrix[0,1],matrix[0].sum()),
        false_exclusion_to_M1=divide(matrix[1,0],matrix[1].sum()),
        false_exclusion_to_UNKNOWN=divide(matrix[1,2],matrix[1].sum()),
        false_exclusion_total=divide(matrix[1,0]+matrix[1,2],matrix[1].sum()),
        UNKNOWN_fraction=float(np.mean(predicted==0)),sample_count=len(actual),
        confusion_matrix=matrix.tolist())
    native_truth,native_prediction=np.asarray(truth),np.asarray(prediction)
    mevo_count=int(np.isin(native_truth,[2,3]).sum());within={}
    for a,b in [(2,3),(3,2)]:
        count=int(np.sum((native_truth==a)&(native_prediction==b)))
        within[f'M{a}_to_M{b}']=dict(count=count,within_mevo_fraction=divide(count,mevo_count),
                                     true_class_fraction=divide(count,int(np.sum(native_truth==a))))
    within['true_mevo_count']=mevo_count
    within['now_correct_binary_count']=sum(within[k]['count'] for k in ['M2_to_M3','M3_to_M2'])
    result['within_mevo']=within
    return result


def evaluate_mevo_boundary(points,truth,prediction):
    return boundary_error(points,collapse_labels_to_mevo(truth,truth=True),
                          collapse_labels_to_mevo(prediction),(1,2))


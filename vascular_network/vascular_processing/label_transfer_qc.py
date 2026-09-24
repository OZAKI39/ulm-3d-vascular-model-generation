"""Held-out semantic scores and explicit production admission gates."""
import numpy as np
from scipy.spatial import cKDTree

GATES=dict(macro_f1=.70,M2_f1=.65,M3_f1=.65,require_baseline_superiority=True,roi_confidence=.60)


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


def pilot_gate(metric,baseline):
    return bool(metric['macro_f1']>=GATES['macro_f1'] and metric['M2_f1']>=GATES['M2_f1'] and metric['M3_f1']>=GATES['M3_f1'] and metric['macro_f1']>baseline['macro_f1'])

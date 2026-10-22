"""Frozen registration ranking and equal multi-donor nearest-neighbor votes."""
import numpy as np
from scipy.spatial import cKDTree
from .similarity_registration import register,transform_points


def select_donors(pool,target,top_k=5):
    ranked=[]
    for donor in pool:
        if donor.geometry.case_id==target.case_id or donor.geometry.side!=target.side:continue
        registration=register(donor.geometry,target)
        ranked.append((donor,registration))
    ranked.sort(key=lambda item:(item[1]['status']!='VALID',-item[1].get('fitness',0),item[1].get('normalized_rmse',1e9),item[0].geometry.case_id))
    selected=[item for item in ranked if item[1]['status']=='VALID'][:top_k]
    for rank,(_,r) in enumerate(selected,1):r['rank']=rank
    return selected,[r for _,r in ranked]


def nearest_baseline(selected,target):
    votes=[]
    for donor,registration in selected:
        points=transform_points(donor.geometry.points,np.asarray(registration['transform']))
        ids=cKDTree(points).query(target.points)[1]
        votes.append(np.eye(3)[donor.labels[ids]-1])
    probabilities=np.mean(votes,axis=0)
    return probabilities.argmax(axis=1)+1


def nearest_support(selected,points):
    """Existing NN queries plus auditable per-donor votes and distance support.

    Binary production votes combine native M2/M3 *before* voting. No target
    labels are accepted, and the caller supplies the frozen QC/ranked donors.
    """
    votes=[];distances=[];normalized=[];identities=[]
    for donor,registration in selected:
        if registration['status']!='VALID':
            raise ValueError('Only registration-QC-valid donors may vote')
        transformed=transform_points(donor.geometry.points,np.asarray(registration['transform']))
        distance,ids=cKDTree(transformed).query(points)
        rmse=float(registration['inlier_rmse'])
        if not np.isfinite(rmse) or rmse<0:raise ValueError('Invalid registration inlier RMSE')
        votes.append(donor.labels[ids]);distances.append(distance)
        normalized.append(distance/max(rmse,1e-8))
        identities.append(donor.geometry.case_id)
    n=len(points)
    native=np.asarray(votes,dtype=np.int8).reshape(len(votes),n)
    binary=np.where(native==3,2,native)
    counts=np.full(n,len(votes),dtype=np.int16)
    p1=np.mean(binary==1,axis=0) if votes else np.zeros(n)
    p2=np.mean(binary==2,axis=0) if votes else np.zeros(n)
    return dict(per_donor_native_vote=native,per_donor_binary_vote=binary,
        per_donor_distance_mm=np.asarray(distances).reshape(len(votes),n),
        per_donor_normalized_distance=np.asarray(normalized).reshape(len(votes),n),
        donor_ids=np.asarray(identities,dtype='U16'),valid_donor_count=counts,
        vote_M1=(binary==1).sum(axis=0),vote_MeVO=(binary==2).sum(axis=0),
        p_M1=p1,p_MeVO=p2,agreement_M1=p1,agreement_MeVO=p2,
        agreement=np.maximum(p1,p2),vote_margin=abs(p2-p1),
        winner=np.where(counts>0,np.where(p1>=p2,1,2),0).astype(np.int8),
        median_normalized_distance=np.median(normalized,axis=0) if votes else np.full(n,np.inf))

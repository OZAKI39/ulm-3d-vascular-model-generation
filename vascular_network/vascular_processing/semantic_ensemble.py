"""Equal donor votes with explicit missing-support accounting."""
import numpy as np
from scipy.spatial import cKDTree
from .similarity_registration import register,transform_points

ENSEMBLE_RULES=dict(top_k=5,weighting='equal among supported donors at each point',minimum_supported_donors_fraction=.5,confidence_gate=.60)


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


def ensemble(results):
    if not results:raise ValueError('No valid donor transports')
    valid=np.array([r['valid'] for r in results]);counts=valid.sum(axis=0)
    probs=np.array([r['probabilities'] for r in results])
    mean=np.divide((probs*valid[:,:,None]).sum(axis=0),counts[:,None],out=np.zeros_like(probs[0]),where=counts[:,None]>0)
    support=np.mean([r['support_mass'] for r in results],axis=0)
    enough=counts>=max(1,int(np.ceil(.5*len(results))))
    confidence=np.where(enough,mean.max(axis=1),0.)
    return dict(probabilities=mean,support_mass=support,valid_donor_count=counts,
                confidence=confidence,predicted_label=np.where(enough,mean.argmax(axis=1)+1,0),valid=enough)


def nearest_baseline(selected,target):
    votes=[]
    for donor,registration in selected:
        points=transform_points(donor.geometry.points,np.asarray(registration['transform']))
        ids=cKDTree(points).query(target.points)[1]
        votes.append(np.eye(3)[donor.labels[ids]-1])
    probabilities=np.mean(votes,axis=0)
    return probabilities.argmax(axis=1)+1

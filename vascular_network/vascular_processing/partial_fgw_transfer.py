"""POT's native Partial Fused Gromov-Wasserstein and soft semantic transport."""
import numpy as np
from scipy.spatial.distance import cdist
from .similarity_registration import transform_points
from .transfer_runtime import require_versions

SUPPORT_RULE = 'column mass >= 0.25 * median(positive column masses); numerical zero excluded'


def normalize_cost(cost):
    values=cost[cost>1e-12]
    scale=float(np.median(values)) if len(values) else 1.
    return cost/scale,scale


def transport_probabilities(transport, labels):
    transport=np.asarray(transport,float)
    if not np.isfinite(transport).all() or transport.min() < -1e-9:raise ValueError('Invalid native transport matrix')
    transport=np.maximum(transport,0.)
    mass=transport.sum(axis=0);positive=mass[mass>1e-12]
    threshold=.25*float(np.median(positive)) if len(positive) else float('inf')
    label_mass=np.array([transport[np.asarray(labels)==c].sum(axis=0) for c in [1,2,3]]).T
    probabilities=np.divide(label_mass,mass[:,None],out=np.zeros_like(label_mass),where=mass[:,None]>1e-12)
    valid=(mass>=threshold)&(mass>1e-12)
    predicted=np.where(valid,probabilities.argmax(axis=1)+1,0)
    return dict(probabilities=probabilities,support_mass=mass,valid=valid,predicted_label=predicted,
                confidence=np.where(valid,probabilities.max(axis=1),0.),support_threshold=threshold,support_rule=SUPPORT_RULE)


def transfer(donor, target, registration, *, alpha=.5, mass_fraction=.8, debug_path=None):
    require_versions()
    import ot
    if registration['status']!='VALID':raise ValueError('Invalid registration cannot enter FGW')
    if donor.geometry.side!=target.side:raise ValueError('Donor side mismatch')
    aligned=transform_points(donor.geometry.points,np.asarray(registration['transform']))
    c1,s1=normalize_cost(cdist(aligned,aligned));c2,s2=normalize_cost(cdist(target.points,target.points))
    feature,sf=normalize_cost(cdist(aligned,target.points))
    p,q=donor.geometry.mass,target.mass;m=mass_fraction*min(p.sum(),q.sum())
    coupling,log=ot.gromov.partial_fused_gromov_wasserstein(feature,c1,c2,p=p,q=q,m=m,alpha=alpha,
                         loss_fun='square_loss',symmetric=True,log=True,numItermax=10000,tol=1e-8)
    if not np.isclose(coupling.sum(),m,atol=1e-6):raise ValueError('POT mass constraint failed')
    if np.any(coupling.sum(axis=1)>p+1e-7) or np.any(coupling.sum(axis=0)>q+1e-7):raise ValueError('POT marginal capacity exceeded')
    result=transport_probabilities(coupling,donor.labels)
    loss=log.get('partial_fgw_dist',log.get('loss',[]))
    if isinstance(loss,(list,np.ndarray)):loss=float(loss[-1]) if len(loss) else None
    elif loss is not None:loss=float(loss)
    result.update(alpha=alpha,mass_fraction=mass_fraction,transported_mass=float(coupling.sum()),loss=loss,
                  iterations=len(log.get('loss',[])),normalization_scales=dict(C1=s1,C2=s2,M=sf),source_case=donor.geometry.case_id)
    if debug_path is not None:np.savez_compressed(debug_path,C1=c1,C2=c2,M=feature,T=coupling)
    return result

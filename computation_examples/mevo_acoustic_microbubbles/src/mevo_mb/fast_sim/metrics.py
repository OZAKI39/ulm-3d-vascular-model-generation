"""Sampled connected components and shape tensors, with fixed particle identities."""
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import squareform,pdist

DEFINITIONS={'cluster_gap_over_radius':0.1,'near_contact_gap_over_radius':0.02,
             'hard_contact_gap_tolerance_over_radius':0.001,
             'meaning':'spatial connected components; no adhesion or fusion'}


def frame_metrics(x,radius):
    n=len(x); distance=squareform(pdist(x)); np.fill_diagonal(distance,np.inf)
    gap=distance-2*radius
    count,labels=connected_components(gap<=DEFINITIONS['cluster_gap_over_radius']*radius,directed=False)
    sizes=np.bincount(labels,minlength=n); centroids=np.zeros((n,3)); tensors=np.zeros((n,3,3))
    eigen=np.zeros((n,3)); axes=np.zeros((n,3,3)); rg=np.zeros(n); identifiable=np.zeros(n,bool)
    for c in range(count):
        xx=x[labels==c]; centroids[c]=xx.mean(axis=0); rel=xx-centroids[c]
        tensors[c]=rel.T@rel/len(xx); e,a=np.linalg.eigh(tensors[c]); eigen[c]=e[::-1]; axes[c]=a[:,::-1]
        rg[c]=np.sqrt(np.trace(tensors[c])); identifiable[c]=bool(e[-1]>0 and e[-1]-e[-2]>1e-3*e[-1])
    largest=int(np.argmax(sizes)); nn=distance.min(axis=1)
    return dict(cluster_count=count,labels=labels,cluster_sizes=sizes,centroids_m=centroids,
                shape_tensors_m2=tensors,eigenvalues_m2=eigen,principal_axes=axes,axis_identifiable=identifiable,
                rg_m=rg,largest_size=int(sizes[largest]),largest_centroid_m=centroids[largest],
                largest_rg_m=rg[largest],largest_eigenvalues_m2=eigen[largest],
                global_centroid_m=x.mean(axis=0),global_rg_m=float(np.sqrt(np.mean(np.sum((x-x.mean(axis=0))**2,axis=1)))),
                mean_nearest_neighbor_m=float(nn.mean()),minimum_gap_m=float(gap.min()),
                near_contact_pair_count=int(np.count_nonzero(np.triu(gap<=0.02*radius,1))),
                hard_contact_pair_count=int(np.count_nonzero(np.triu(gap<=0.001*radius,1))))


def trajectory_metrics(positions,radius,times):
    records=[[frame_metrics(x,radius) for x in frame] for frame in positions]
    arrays={key:np.asarray([[r[key] for r in frame] for frame in records]) for key in records[0][0]}
    events=[]
    for t in range(1,len(times)):
        for b in range(positions.shape[1]):
            old=arrays['labels'][t-1,b]; new=arrays['labels'][t,b]
            for kind,left,right in [('merge',old,new),('split',new,old)]:
                for c in np.unique(right):
                    members=np.flatnonzero(right==c); others=np.unique(left[members])
                    if len(others)>1:
                        events.append(dict(kind=kind,batch=b,observed_time_s=float(times[t]),
                            bracket_s=[float(times[t-1]),float(times[t])],particle_ids=members.tolist(),
                            scope='sampled membership change; event time bounded by saved frames'))
    return arrays,events

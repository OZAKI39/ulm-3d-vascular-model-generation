"""Frozen engineering categories; uncovered cases remain explicitly uncertain."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

def clot_components(n,bonds,types,alive):
    e=bonds[alive&(types==1)]
    graph=coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(n,n)).tocsr()
    nc,lab=connected_components(graph);counts=np.bincount(lab,minlength=nc)
    attached=np.zeros(nc,bool);attached[lab[bonds[alive&(types==2),0]]]=True
    return lab,counts,attached

def classify(m):
    if m['debris_fraction']>=.5 or (m['median_fragment_nodes']<=2 and m['bulk_broken_fraction']>=.5):return 'PULVERIZATION'
    if m['interface_broken_fraction']>=.9 and m['largest_free_fraction']>=.8 and (m['resolved_fragment_count']==1 or m['bulk_broken_fraction']<.01):return 'WHOLE_CLOT_DETACHMENT'
    if m['resolved_fragment_count']<=1 and m['interface_broken_fraction']<.1 and m['bulk_broken_fraction']<.01:return 'NO_FAILURE'
    if m['interface_broken_fraction']>=.5 and m['resolved_fragment_count']>=2:return 'HYBRID_DEBONDING_AND_FRAGMENTATION'
    if m['resolved_fragment_count']>=2 and m['debris_fraction']<.5:return 'RESOLVED_FRAGMENTATION'
    return 'CLASSIFICATION_UNCERTAIN'

"""Only surviving tensile bonds, separately for bulk and interface."""
import numpy as np

def utilization(X,bonds,types,alive,L0,sbulk,sinterface):
    out={}
    for typ,name,critical in [(1,'bulk',sbulk),(2,'interface',sinterface)]:
        mask=alive&(types==typ);e=bonds[mask]
        u=np.maximum(np.linalg.norm(X[e[:,1]]-X[e[:,0]],axis=1)/L0[mask]-1,0)/critical if len(e) else np.empty(0)
        out['surviving_'+name+'_bonds']=int(len(u))
        for q,key in [(99,'P99'),(99.9,'P99_9')]:out['U_'+name+'_'+key]=float(np.percentile(u,q)) if len(u) else None
        for threshold,key in [(.8,'0_8'),(.9,'0_9')]:out['frac_'+name+'_U_gt_'+key]=float(np.mean(u>threshold)) if len(u) else None
    return out

"""Physical sample locations, volume weights and frozen ROI; no field fitting."""
import numpy as np
from scipy.spatial import cKDTree
from .geometry import ring_distance

def distribution(values,weights,points,mask=None):
    mask=np.ones(len(values),bool) if mask is None else mask
    ids=np.flatnonzero(mask);v=np.asarray(values)[ids];w=np.asarray(weights)[ids]
    if not len(ids):raise ValueError('Empty declared ROI')
    order=np.argsort(v,kind='stable');c=np.cumsum(w[order]);result={}
    for name,q in [('P95',.95),('P99',.99)]:
        i=ids[order[min(np.searchsorted(c,q*c[-1],side='left'),len(order)-1)]]
        result[name]=float(values[i]);result[name+'_index']=int(i);result[name+'_position_m']=points[i]
    i=ids[np.argmax(v)];result.update(max=float(values[i]),max_index=int(i),max_position_m=points[i],weighted_RMS=float(np.sqrt(np.average(v*v,weights=w))),samples=len(v),weight_sum=float(w.sum()))
    return result

def summary_fields(X,V,pairs,bond_weight,u,F,Ft,psin,psib,config):
    mid=(X[pairs[:,0]]+X[pairs[:,1]])/2;R=config['R_m'];width=config['sampling']['ROI_width_m']
    node_roi=ring_distance(X,R)>=width;bond_roi=ring_distance(mid,R)>=width
    apex=int(np.argmin(np.linalg.norm(X-config['sampling']['apex'],axis=1)))
    result=dict(apex_displacement_m=u[apex],apex_position_m=X[apex],u_RMS_m=np.sqrt(np.average(np.sum(u*u,axis=1),weights=V)))
    for name,value,weight,points,roi in [('node_Psi_plus',psin,V,X,node_roi),('bond_Psi_plus',psib,bond_weight,mid,bond_roi),('Fbar_minus_I',np.linalg.norm(F-np.eye(3),axis=(1,2)),V,X,node_roi),('Ftilde_minus_I',np.linalg.norm(Ft-np.eye(3),axis=(1,2)),bond_weight,mid,bond_roi)]:
        result[name]=dict(full=distribution(value,weight,points),ROI=distribution(value,weight,points,roi))
    return result

def profile_points(config):
    a=config['sampling']['profiles']['axis'];b=config['sampling']['profiles']['transverse']
    axis=np.zeros((a['n'],3));axis[:,2]=np.linspace(a['z_min_m'],a['z_max_m'],a['n'])
    cross=np.zeros((b['n'],3));cross[:,0]=np.linspace(b['x_min_m'],b['x_max_m'],b['n']);cross[:,2]=b['z_m']
    return dict(axis=axis,transverse=cross)

def interpolate_tets(X,cells,values,points):
    """P1 interpolation using the actual supplied cells, never a new triangulation."""
    centers=X[cells].mean(1);tree=cKDTree(centers);result=[]
    for q in points:
        found=False
        for k in (min(64,len(cells)),len(cells)):
            ids=np.atleast_1d(tree.query(q,k=k)[1]);cc=cells[ids]
            edge=(X[cc[:,1:]]-X[cc[:,0,None]]).swapaxes(1,2)
            bary=np.linalg.solve(edge,(q-X[cc[:,0]])[...,None])[...,0];bary=np.column_stack([1-bary.sum(1),bary])
            hit=np.flatnonzero(np.min(bary,axis=1)>=-1e-9)
            if len(hit):
                j=hit[0];result.append(np.tensordot(bary[j],values[cc[j]],axes=(0,0)));found=True;break
        if not found:raise ValueError(f'Profile point outside saved actual mesh: {q}')
    return np.asarray(result)

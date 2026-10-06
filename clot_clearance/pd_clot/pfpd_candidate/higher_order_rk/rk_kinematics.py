"""Displacement output index first; physical derivative index second."""
import numpy as np

def gradient(values,source,target,g):
    values=np.asarray(values)
    if values.ndim==1:values=values[:,None]
    du=values[target]-values[source];out=np.empty((len(values),values.shape[1],3))
    for a in range(values.shape[1]):
        for c in range(3):out[:,a,c]=np.bincount(source,weights=du[:,a]*g[:,c],minlength=len(values))
    return out

def deformation(X,x,source,target,g):return np.eye(3)+gradient(x-X,source,target,g)

def c1_baseline(X,x,V,pairs,omega):
    """Unchanged baqp.evaluate nodal algebra only, excluding stress, Z and force.

    K solve replaces explicit inverse algebraically; regression checks original saved F.
    """
    i,j=pairs.T;xi=X[j]-X[i];n=len(X)
    K=np.zeros((n,3,3));q=omega[:,None,None]*xi[:,:,None]*xi[:,None,:]
    np.add.at(K,i,q*V[j,None,None]);np.add.at(K,j,q*V[i,None,None])
    solved=np.linalg.solve(K,np.broadcast_to(np.eye(3),K.shape))
    Bi=omega[:,None]*np.einsum('bij,bj->bi',solved[i],xi)
    Bj=-omega[:,None]*np.einsum('bij,bj->bi',solved[j],xi)
    du=(x-X)[j]-(x-X)[i];F=np.broadcast_to(np.eye(3),(n,3,3)).copy()
    np.add.at(F,i,du[:,:,None]*Bi[:,None,:]*V[j,None,None])
    np.add.at(F,j,-du[:,:,None]*Bj[:,None,:]*V[i,None,None])
    return F

"""Kinematic midpoint projection and unchanged positive-energy density only."""
import numpy as np
from ..svk_spectral_split import green_lagrange,svk_split_energy

def associated_gradient(X,x,pairs,F):
    i,j=pairs.T;xi=X[j]-X[i];eta=x[j]-x[i];av=.5*(F[i]+F[j])
    residual=eta-np.einsum('bij,bj->bi',av,xi)
    return av+residual[:,:,None]*xi[:,None,:]/np.sum(xi*xi,axis=1)[:,None,None]

def psi(F):return svk_split_energy(green_lagrange(F),0.,1000.)[0]

def metrics(error,weight,scale):
    e=np.asarray(error);value=np.sqrt(np.sum(e*e,axis=tuple(range(1,e.ndim)))) if e.ndim>1 else np.abs(e)
    return dict(max=float(value.max()),weighted_RMS=float(np.sqrt(np.average(value*value,weights=weight))),
        max_scaled=float(value.max()/scale),weighted_RMS_scaled=float(np.sqrt(np.average(value*value,weights=weight))/scale))

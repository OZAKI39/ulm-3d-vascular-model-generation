"""BAQP fixed-state mechanics, Eq2.10–2.19 with author-code r² correction.

No time integration, no phase evolution, no singular-family regularization.
"""
from dataclasses import dataclass
import numpy as np
from .energy import Material,constitutive


@dataclass
class State:
    K: np.ndarray
    Fbar: np.ndarray
    Ftilde: np.ndarray
    shape_i: np.ndarray
    shape_j: np.ndarray
    W0: np.ndarray
    P: np.ndarray
    Z: np.ndarray
    point_energy: np.ndarray
    bond_energy: np.ndarray
    force: np.ndarray
    force_state_i: np.ndarray
    force_state_j: np.ndarray

    @property
    def energy(self):return float(self.bond_energy.sum())


def evaluate(family,x,material=Material(),*,g=None,h=None):
    f=family;x=np.asarray(x,float);i,j=f.pairs.T;V=f.volume;m=len(i);n=len(V)
    if x.shape!=f.X.shape or not np.isfinite(x).all():raise ValueError('Invalid coordinates')
    def factor(v):
        a=np.ones(m) if v is None else np.asarray(v,float)
        if a.shape!=(m,) or np.any((a<0)|(a>1)) or not np.isfinite(a).all():raise ValueError('Invalid prescribed factor')
        return a
    g=factor(g);h=factor(h);w=f.omega*h;xi=f.xi
    K=np.zeros((n,3,3));q=w[:,None,None]*xi[:,:,None]*xi[:,None,:]
    np.add.at(K,i,q*V[j,None,None]);np.add.at(K,j,q*V[i,None,None])
    eig=np.linalg.eigvalsh(K)
    if np.any(eig[:,0]<=np.maximum(eig[:,-1],0)*1e-12):raise ValueError('Rank deficient family; outside candidate scope')
    inv=np.linalg.inv(K)
    Bi=w[:,None]*np.einsum('bij,bj->bi',inv[i],xi)
    Bj=-w[:,None]*np.einsum('bij,bj->bi',inv[j],xi)
    # Displacement form as official rkc_defgrad!: affine completeness gives
    # exact equivalence to the paper's total-position form for full-rank K.
    du=(x-f.X)[j]-(x-f.X)[i]
    F=np.broadcast_to(np.eye(3),(n,3,3)).copy()
    np.add.at(F,i,du[:,:,None]*Bi[:,None,:]*V[j,None,None])
    np.add.at(F,j,-du[:,:,None]*Bj[:,None,:]*V[i,None,None])
    av=.5*(F[i]+F[j]);eta=x[j]-x[i];r2=f.length**2
    residual=eta-np.einsum('bij,bj->bi',av,xi)
    Ft=av+residual[:,:,None]*xi[:,None,:]/r2[:,None,None]
    W,P0=constitutive(Ft,material);P=g[:,None,None]*P0
    factor_pair=f.wbar*V[i]*V[j]
    bond_energy=2*factor_pair*g*W
    point_energy=np.zeros(n)
    np.add.at(point_energy,i,factor_pair*g*W);np.add.at(point_energy,j,factor_pair*g*W)
    proj=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r2[:,None,None]
    Z=np.zeros((n,3,3));term=f.wbar[:,None,None]*(P@proj)
    np.add.at(Z,i,term*V[j,None,None]);np.add.at(Z,j,term*V[i,None,None])
    direct=f.wbar[:,None]*np.einsum('bij,bj->bi',P,xi)/r2[:,None]
    ti=direct+np.einsum('bij,bj->bi',Z[i],Bi)
    tj=-direct+np.einsum('bij,bj->bi',Z[j],Bj)
    pairforce=(ti-tj)*V[i,None]*V[j,None]
    force=np.zeros((n,3));np.add.at(force,i,pairforce);np.add.at(force,j,-pairforce)
    return State(K,F,Ft,Bi,Bj,W,P,Z,point_energy,bond_energy,force,ti,tj)

"""The unchanged intact BAQP bond potential, with one physical-bond prefactor."""
import numpy as np
from ..energy import constitutive


def bond_prefactors(family):
    i,j=family.pairs.T
    return 2*(family.wbar*family.volume[i]*family.volume[j])


def bond_state(reference,x,F,sl,material):
    i,j=reference.pairs[sl].T
    xi=reference.xi[sl];r2=reference.r2[sl]
    avg=.5*(F[i]+F[j]);eta=x[j]-x[i]
    mismatch=eta-np.einsum('bij,bj->bi',avg,xi)
    Ft=avg+mismatch[:,:,None]*xi[:,None,:]/r2[:,None,None]
    W,P=constitutive(Ft,material)
    return Ft,W,P


def energy_from_nodes(reference,x,F,material,block=20000):
    total=0.
    for start in range(0,len(reference.pairs),block):
        sl=slice(start,start+block)
        _,W,_=bond_state(reference,x,F,sl,material)
        total+=float(np.sum(reference.m[sl]*W))
    return total

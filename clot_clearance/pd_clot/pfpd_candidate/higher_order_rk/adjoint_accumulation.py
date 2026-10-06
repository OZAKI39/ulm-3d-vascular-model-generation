"""Linear-work two-stage adjoint; no per-bond neighbor-of-neighbor loops."""
import numpy as np


def bond_adjoint(reference,P,sl,G,gradU):
    i,j=reference.pairs[sl].T;xi=reference.xi[sl];r2=reference.r2[sl]
    R=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r2[:,None,None]
    A=.5*reference.m[sl,None,None]*(P@R)
    t=reference.m[sl,None]*np.einsum('bij,bj->bi',P,xi)/r2[:,None]
    # Direct indexed accumulation avoids a length-N buffer for every bond block.
    np.add.at(G,i,A);np.add.at(G,j,A)
    np.add.at(gradU,j,t);np.add.at(gradU,i,-t)


def nodal_adjoint(operator,G,gradU,block=100000):
    for start in range(0,len(operator.source),block):
        sl=slice(start,start+block);i=operator.source[sl];k=operator.target[sl]
        q=np.einsum('bij,bj->bi',G[i],operator.g[sl])
        np.add.at(gradU,k,q);np.add.at(gradU,i,-q)
    return gradU


def delta_bond_gradient(reference,delta_x,delta_F,sl):
    i,j=reference.pairs[sl].T;xi=reference.xi[sl];r2=reference.r2[sl]
    R=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r2[:,None,None]
    return .5*(delta_F[i]+delta_F[j])@R+(delta_x[j]-delta_x[i])[:,:,None]*xi[:,None,:]/r2[:,None,None]

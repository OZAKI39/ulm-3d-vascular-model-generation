"""Constitutive energy and PK stress at bond quadrature points; SI units."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Material:
    model: str='svk'
    mu: float=1000.
    bulk: float=50000.

    def __post_init__(self):
        if self.model not in ('svk','neo_hookean') or not np.isfinite([self.mu,self.bulk]).all() or min(self.mu,self.bulk)<=0:raise ValueError('Invalid constitutive parameters')

    @property
    def lame(self):return self.bulk-2*self.mu/3


def constitutive(F,material):
    F=np.asarray(F,float);I=np.eye(3);mu=material.mu;lam=material.lame
    J=np.linalg.det(F)
    if not np.isfinite(F).all() or np.any(J<=0):raise ValueError('Positive determinant required')
    if material.model=='svk':
        E=.5*(np.swapaxes(F,-1,-2)@F-I);tr=np.trace(E,axis1=-2,axis2=-1)
        W=.5*lam*tr**2+mu*np.sum(E*E,axis=(-1,-2))
        P=F@(lam*tr[...,None,None]*I+2*mu*E)
    else:
        # Same compressible NH potential as production. Singular values + log1p
        # avoid subtracting F:F-3 at tiny historical strains; no energy clipping.
        d=np.linalg.svd(F,compute_uv=False)-1
        l=np.log1p(d);logJ=l.sum(axis=-1)
        W=.5*mu*np.sum(d*d+2*(d-l),axis=-1)+.5*lam*logJ**2
        invT=np.swapaxes(np.linalg.inv(F),-1,-2)
        P=mu*(F-invT)+lam*logJ[...,None,None]*invT
    return W,P

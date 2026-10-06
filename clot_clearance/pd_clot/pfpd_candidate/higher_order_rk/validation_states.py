"""Deterministic validation fields; reuses the accepted manufactured targets."""
import numpy as np
from .manufactured_fields import exact


def rotation(angle_degrees=37.,axis=(1.,2.,3.)):
    a=np.array(axis,float);a/=np.linalg.norm(a)
    K=np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0.]])
    t=np.deg2rad(angle_degrees)
    return np.eye(3)+np.sin(t)*K+(1-np.cos(t))*(K@K)


def smooth_direction(X,L,seed):
    q=X/L;coeff=np.random.default_rng(seed).normal(size=(9,3))
    functions=np.column_stack([q,np.sin(np.pi*q),np.cos(np.pi*q)])
    v=functions@coeff;v-=v.mean(axis=0);v/=np.linalg.norm(v)
    return v


def target(field,X,L):
    if field=='finite_affine_shear':
        F=np.eye(3);F[0,2]=.03
        return X@F.T,np.broadcast_to(F,(len(X),3,3))
    return exact(field,X,L)

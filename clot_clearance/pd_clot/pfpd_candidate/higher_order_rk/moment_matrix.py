"""Reference local moment systems. No inverse, pseudo-inverse or regularization."""
import numpy as np
from .basis import basis,selector

def local_system(xi,omega,neighbor_volume,delta,dimension=9):
    Q=basis(xi/delta,dimension);w=omega*neighbor_volume
    M=Q.T@(w[:,None]*Q)
    singular=np.linalg.svd(M,compute_uv=False)
    tolerance=dimension*np.finfo(float).eps*singular[0]
    rank=int(np.count_nonzero(singular>tolerance))
    condition=float(singular[0]/singular[-1]) if singular[-1]>0 else np.inf
    coeff=np.full((dimension,3),np.nan)
    if rank==dimension:coeff=np.linalg.solve(M,selector(delta,dimension))
    g=w[:,None]*(Q@coeff)
    return M,singular,rank,condition,coeff,g

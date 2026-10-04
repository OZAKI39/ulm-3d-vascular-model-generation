"""Scalar-loop copy of the existing full-rank NOSB kernel, same energy/forces.

Avoids per-bond small BLAS allocations; reduced-rank cases use the original
explicit intrinsic continuation. No constitutive or support rule is changed.
"""
import numpy as np
from numba import njit
from ..fragment_mechanics import evaluate_supported as reference


@njit(cache=True)
def kernel(x,pairs,xi,weight,g,volume,inv,trace,remaining,shear,bulk,alpha,minimum_J):
    n=len(x);A=np.zeros((n,3,3));F=np.zeros_like(A);PK=np.zeros_like(A)
    E=np.zeros(n);H=np.zeros(n);J=np.zeros(n);coefficient=np.zeros(n)
    for k in range(len(pairs)):
        i,j=pairs[k];q=weight[k]*g[k]
        for a in range(3):
            eta=x[j,a]-x[i,a]
            for b in range(3):
                v=q*eta*xi[k,b]
                A[i,a,b]+=v*volume[j];A[j,a,b]+=v*volume[i]
    lame=bulk-2*shear/3
    for i in range(n):
        for a in range(3):
            for b in range(3):
                for d in range(3):F[i,a,b]+=A[i,a,d]*inv[i,d,b]
        f=F[i];cof=np.empty((3,3))
        cof[0,0]=f[1,1]*f[2,2]-f[1,2]*f[2,1]
        cof[0,1]=f[1,2]*f[2,0]-f[1,0]*f[2,2]
        cof[0,2]=f[1,0]*f[2,1]-f[1,1]*f[2,0]
        cof[1,0]=f[0,2]*f[2,1]-f[0,1]*f[2,2]
        cof[1,1]=f[0,0]*f[2,2]-f[0,2]*f[2,0]
        cof[1,2]=f[0,1]*f[2,0]-f[0,0]*f[2,1]
        cof[2,0]=f[0,1]*f[1,2]-f[0,2]*f[1,1]
        cof[2,1]=f[0,2]*f[1,0]-f[0,0]*f[1,2]
        cof[2,2]=f[0,0]*f[1,1]-f[0,1]*f[1,0]
        jac=f[0,0]*cof[0,0]+f[0,1]*cof[0,1]+f[0,2]*cof[0,2];J[i]=jac
        if not np.isfinite(jac) or jac<=minimum_J:raise ValueError('Invalid full-rank deformation Jacobian')
        logJ=np.log(jac);squared=0.
        for a in range(3):
            for b in range(3):
                squared+=f[a,b]**2
                for d in range(3):
                    P=shear*f[a,d]+(lame*logJ-shear)*cof[a,d]/jac
                    PK[i,a,b]+=remaining[i]*P*inv[i,d,b]
        E[i]=remaining[i]*(.5*shear*(squared-3)-shear*logJ+.5*lame*logJ**2)
        coefficient[i]=alpha*shear*remaining[i]/trace[i]
    force=np.zeros((n,3))
    for k in range(len(pairs)):
        q=weight[k]*g[k]
        if q==0:continue
        i,j=pairs[k];ri2=0.;rj2=0.
        for a in range(3):
            ri=x[j,a]-x[i,a];rj=ri;elastic=0.
            for b in range(3):
                ri-=F[i,a,b]*xi[k,b];rj-=F[j,a,b]*xi[k,b]
                elastic+=(PK[i,a,b]+PK[j,a,b])*xi[k,b]
            fij=q*volume[i]*volume[j]*(elastic+coefficient[i]*ri+coefficient[j]*rj)
            force[i,a]+=fij;force[j,a]-=fij
            ri2+=ri*ri;rj2+=rj*rj
        H[i]+=.5*coefficient[i]*q*ri2*volume[j]
        H[j]+=.5*coefficient[j]*q*rj2*volume[i]
    return force,F,J,E,H


def evaluate_supported(cloud,x,integrity,material,safety,shape):
    if not np.all(shape[3]==3):
        return reference(cloud,x,integrity,material,safety,shape)
    return kernel(x,cloud.pairs,cloud.xi,cloud.weight,integrity,cloud.volume,*shape[:3],
                  material['shear_modulus_Pa'],material['bulk_modulus_Pa'],material['stabilization_alpha'],safety['minimum_J'])

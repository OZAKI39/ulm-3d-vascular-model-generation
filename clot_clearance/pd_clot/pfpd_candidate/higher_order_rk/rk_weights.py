"""Cache physical directed gradient weights including V_j, units 1/m."""
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from .moment_matrix import local_system

def directed_topology(pairs,n):
    i,j=pairs.T;source=np.r_[i,j];target=np.r_[j,i]
    order=np.lexsort((target,source));source=source[order];target=target[order]
    ptr=np.r_[0,np.cumsum(np.bincount(source,minlength=n))]
    return source,target,ptr

def build_weights(X,V,source,target,ptr,omega,delta,dimension=9,workers=8,block=128):
    n=len(X);g=np.empty((len(source),3));M=np.empty((n,dimension,dimension))
    sv=np.empty((n,dimension));rank=np.empty(n,dtype=np.int16);cond=np.empty(n);coeff=np.empty((n,dimension,3))
    def work(bounds):
        lo,hi=bounds
        for i in range(lo,hi):
            sl=slice(ptr[i],ptr[i+1]);j=target[sl]
            a,b,c,d,e,f=local_system(X[j]-X[i],omega[sl],V[j],delta,dimension)
            M[i]=a;sv[i]=b;rank[i]=c;cond[i]=d;coeff[i]=e;g[sl]=f
    with ThreadPoolExecutor(workers) as pool:list(pool.map(work,[(k,min(n,k+block)) for k in range(0,n,block)]))
    return dict(M=M,singular_values=sv,rank=rank,condition=cond,coefficients=coeff,g=g)

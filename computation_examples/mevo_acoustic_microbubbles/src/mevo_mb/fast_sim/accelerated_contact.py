"""FP64 mathematical equivalents: all pairs retained, only representation changes.

B has six nonzeros per row. Products evaluate n.(v_j-v_i) and scatter
(-lambda*n,+lambda*n). No distance filtering of contact inequalities occurs.
Resistance uses exactly the historical activation, floor and phi formula.
"""
from functools import lru_cache
import numpy as np
from scipy import linalg

@lru_cache(maxsize=16)
def indices(n):
    i,j=np.triu_indices(n,1)
    i.flags.writeable=False;j.flags.writeable=False
    return i,j

def normals(x):
    x=np.asarray(x,dtype=float)
    i,j=indices(len(x));dx=x[j]-x[i];d=np.linalg.norm(dx,axis=1)
    if np.any(d<=0):raise ValueError('Coincident centers are not a valid contact normal')
    return i,j,d,dx/d[:,None]

class Transpose:
    def __init__(self,base):self.base=base
    def __matmul__(self,lam):
        b=self.base;f=np.asarray(lam)[:,None]*b.normals
        out=np.zeros((b.n,3));np.add.at(out,b.i,-f);np.add.at(out,b.j,f)
        return out.ravel()

class PairMatrix:
    def __init__(self,n,i,j,normals):self.n=n;self.i=i;self.j=j;self.normals=normals;self.shape=(len(i),3*n)
    def __len__(self):return len(self.i)
    @property
    def T(self):return Transpose(self)
    def __matmul__(self,v):
        q=np.asarray(v).reshape(self.n,3)
        return np.einsum('ij,ij->i',self.normals,q[self.j]-q[self.i])
    def rows(self,ids):
        ids=np.asarray(ids,dtype=int);out=np.zeros((len(ids),3*self.n));r=np.arange(len(ids))[:,None];k=np.arange(3)
        out[r,3*self.i[ids,None]+k]=-self.normals[ids]
        out[r,3*self.j[ids,None]+k]=self.normals[ids]
        return out
    def __getitem__(self,ids):return self.rows(ids)
    @property
    def nbytes(self):return self.i.nbytes+self.j.nbytes+self.normals.nbytes

def contact_geometry(x,radius):
    i,j,d,n=normals(x)
    return PairMatrix(len(x),i,j,n),np.column_stack((i,j)),d-2*radius

def resistance(x,radius,activation=.1,floor=.001):
    i,j,d,n=normals(x);e=(d-2*radius)/radius;mask=e<activation
    i,j,n,e=i[mask],j[mask],n[mask],e[mask]
    def phi(v):return .5/v-.45*np.log(v)-(3/56)*v*np.log(v)
    c=.5*(phi(np.maximum(e,floor))-phi(activation))
    A=np.eye(3*len(x));blocks=c[:,None,None]*n[:,:,None]*n[:,None,:]
    axis=np.arange(3);u=3*i[:,None]+axis;v=3*j[:,None]+axis
    # Off-diagonal blocks are unique; diagonal repeated indices must accumulate.
    np.add.at(A,(u[:,:,None],u[:,None,:]),blocks)
    np.add.at(A,(v[:,:,None],v[:,None,:]),blocks)
    A[u[:,:,None],v[:,None,:]]-=blocks;A[v[:,:,None],u[:,None,:]]-=blocks
    return dict(A=A,pairs=np.column_stack((i,j)),normals=n,coefficients=c)

class CPUFactor:
    def __init__(self,A):self.factor=linalg.cho_factor(A,lower=True);self.solve_columns=0
    def solve(self,b):
        self.solve_columns+=1 if np.ndim(b)==1 else b.shape[1]
        return linalg.cho_solve(self.factor,b)

class GPUFactor:
    """Deliberate hybrid comparison, transfer-inclusive; not claimed device resident."""
    def __init__(self,A):
        import cupy as cp
        import cupyx
        self.cp=cp;self.solve_columns=0
        with cupyx.errstate(linalg='raise'):self.L=cp.linalg.cholesky(cp.asarray(A,dtype=cp.float64))
        cp.cuda.get_current_stream().synchronize()
    def solve(self,b):
        from cupyx.scipy.linalg import solve_triangular
        cp=self.cp;self.solve_columns+=1 if np.ndim(b)==1 else b.shape[1]
        rhs=cp.asarray(b,dtype=cp.float64)
        y=solve_triangular(self.L,rhs,lower=True)
        return solve_triangular(self.L.T,y,lower=False).get()

def factorize(A,backend='cpu'):
    if backend=='cpu':return CPUFactor(A)
    if backend=='cupy':return GPUFactor(A)
    raise ValueError(backend)

def active_set_qp(A,factor,b,B,gap_over_R,max_iterations=96):
    """Original working-set decisions/tolerances; lazily solve active columns.

    The cache lives only in THIS QP: A and current B are fixed. It never crosses
    a physical step or an SQP relinearization. Every inactive inequality is checked.
    """
    free=factor.solve(b);active=set(np.flatnonzero(B@free+gap_over_R<0).tolist())
    visited=set();lam=np.zeros(len(B));cache={}
    for iteration in range(max_iterations):
        key=tuple(sorted(active))
        if key in visited:return free,lam,'active_set_cycle',iteration
        visited.add(key);ids=np.array(key,dtype=int);lam=np.zeros(len(B))
        if len(ids):
            missing=[int(k) for k in ids if int(k) not in cache]
            if missing:
                solved=factor.solve(B.rows(missing).T)
                for j,k in enumerate(missing):cache[k]=solved[:,j]
            inverse_BT=np.column_stack([cache[int(k)] for k in ids])
            rows=B.rows(ids);S=rows@inverse_BT;rhs=-gap_over_R[ids]-rows@free
            lam[ids]=linalg.lstsq(S,rhs,cond=1e-12,lapack_driver='gelsy')[0]
            s=free+inverse_BT@lam[ids]
        else:s=free.copy()
        slack=B@s+gap_over_R
        if lam.min(initial=0)<-1e-12:active.remove(int(np.argmin(lam)));continue
        inactive=np.ones(len(B),dtype=bool);inactive[ids]=False
        violations=np.where(inactive,slack,np.inf)
        if violations.min(initial=0)<-1e-12:active.add(int(np.argmin(violations)));continue
        if abs(slack[ids]).max(initial=0)>1e-10:return s,lam,'active_set_inconsistent_equalities',iteration+1
        return s,lam,'active_set_solved',iteration+1
    return s,lam,'active_set_iteration_limit',max_iterations

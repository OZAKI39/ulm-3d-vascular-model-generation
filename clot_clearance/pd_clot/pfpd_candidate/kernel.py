"""One unordered physical bond, two endpoint family index entries (no copied state)."""
from dataclasses import dataclass
import numpy as np


def radial(rho, kind='quadratic'):
    r=np.asarray(rho,dtype=float)
    if np.any(r<0) or not np.isfinite(r).all():raise ValueError('Invalid radius')
    if kind not in ('constant','linear','quadratic'):raise ValueError('Unknown kernel')
    power={'constant':0,'linear':1,'quadratic':2}[kind]
    return np.where(r<1,(np.maximum(1-r,0))**power,0.)


def readonly(a,dtype=float):
    a=np.array(a,dtype=dtype,copy=True);a.setflags(write=False);return a


@dataclass(frozen=True)
class Family:
    X: np.ndarray
    volume: np.ndarray
    pairs: np.ndarray
    xi: np.ndarray
    length: np.ndarray
    omega: np.ndarray
    omega0: np.ndarray
    wbar: np.ndarray
    horizon: float
    # These are indices into the unique physical bonds, not two s arrays.
    oriented_bond: np.ndarray
    source: np.ndarray
    target: np.ndarray
    sign: np.ndarray

    @classmethod
    def from_cloud(cls,cloud):
        return cls.from_arrays(cloud.X,cloud.volume,cloud.pairs,cloud.horizon,cloud.weight)

    @classmethod
    def from_arrays(cls,X,volume,pairs,horizon,omega=None):
        X=readonly(X);V=readonly(volume);p=readonly(pairs,int)
        if X.ndim!=2 or X.shape[1]!=3 or V.shape!=(len(X),) or np.any(V<=0) or not np.isfinite(X).all() or not np.isfinite(V).all():raise ValueError('Invalid points or volume')
        if p.ndim!=2 or p.shape[1]!=2 or np.any(p<0) or np.any(p>=len(X)) or np.any(p[:,0]>=p[:,1]) or len(np.unique(p,axis=0))!=len(p):raise ValueError('Unique i<j physical bonds required')
        if not np.isfinite(horizon) or horizon<=0:raise ValueError('Invalid horizon')
        i,j=p.T;xi=X[j]-X[i];r=np.linalg.norm(xi,axis=1)
        if np.any(r<=0) or np.any(r>horizon*(1+1e-12)):raise ValueError('Invalid bond length')
        w=radial(r/horizon) if omega is None else np.array(omega,float,copy=True)
        if w.shape!=r.shape or np.any(w<0) or not np.isfinite(w).all():raise ValueError('Invalid influence')
        w0=np.zeros(len(X));np.add.at(w0,i,w*V[j]);np.add.at(w0,j,w*V[i])
        if np.any(w0<=0):raise ValueError('Empty weighted family')
        wb=.5*w*(1/w0[i]+1/w0[j]);m=len(p)
        return cls(X,V,p,readonly(xi),readonly(r),readonly(w),readonly(w0),readonly(wb),float(horizon),readonly(np.tile(np.arange(m),2),int),readonly(np.r_[i,j],int),readonly(np.r_[j,i],int),readonly(np.r_[np.ones(m),-np.ones(m)]))

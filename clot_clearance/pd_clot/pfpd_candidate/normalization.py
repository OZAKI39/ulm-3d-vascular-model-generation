"""Continuum c0 and discrete planar crack quadrature. No fracture evolution."""
import numpy as np
from scipy.integrate import quad
from scipy.signal import fftconvolve
from .kernel import radial


def c0(kind='quadratic'):
    numerator=quad(lambda r:float(radial(r,kind))*r**3,0,1,epsabs=1e-13)[0]
    denominator=2*quad(lambda r:float(radial(r,kind))*r**2,0,1,epsabs=1e-13)[0]
    return dict(kernel=kind,numerator=numerator,denominator=denominator,c0=numerator/denominator,analytic={'constant':3/8,'linear':3/10,'quadratic':1/4}[kind])


def crack_energy(family,s,Gc=1.,normalizer=.25):
    from .degradation import bounded
    s=bounded(s)
    if s.shape!=(len(family.pairs),) or Gc<=0 or normalizer<=0:raise ValueError('Invalid crack energy input')
    i,j=family.pairs.T
    return float(np.sum(2*family.wbar*family.volume[i]*family.volume[j]*Gc*s/(normalizer*family.horizon)))


def lattice_stencil(ratio):
    if ratio<=1:raise ValueError('Ratio must exceed 1')
    m=int(np.ceil(ratio));offs=np.indices((2*m+1,)*3).reshape(3,-1).T-m
    r=np.linalg.norm(offs,axis=1);valid=(r>0)&(r<ratio)
    o=offs[valid];w=radial(r[valid]/ratio)
    return o,w,m


def bulk_planar_ratio(ratio):
    """Exact infinite regular-lattice pair counting per unit crack area.

    For offset oz>0, there are oz crossing start layers. Directed double
    sum contributes 2; exclude self from omega0 as in the point cloud.
    """
    o,w,_=lattice_stencil(ratio);z=o[:,2]>0
    return float(2*np.sum(w[z]*o[z,2])/(w.sum()*.25*ratio))


def coupon(n,ratio,side=.003,Gc=1.):
    """Finite n³ regular coupon, plane midway between the two center layers.

    Fast offset enumeration is algebraically the same unordered pair sum.
    omega0 includes the intact truncated families even across the crack.
    No fitted surface correction. Boundary bias is explicitly retained.
    """
    if n%2 or n<4:raise ValueError('Even n required')
    h=side/n;V=h**3;delta=ratio*h;o,w,m=lattice_stencil(ratio)
    stencil=np.zeros((2*m+1,)*3);stencil[tuple((o+m).T)]=w
    omega0=fftconvolve(np.ones((n,n,n)),stencil,mode='same')*V
    accum=0.;crossing=0
    for d,weight in zip(o,w):
        dx,dy,dz=map(int,d)
        if dz<=0 or dz>=n:continue
        sx=slice(max(0,-dx),min(n,n-dx));sy=slice(max(0,-dy),min(n,n-dy))
        z0=max(0,n//2-dz);z1=min(n//2,n-dz)
        if z1<=z0:continue
        a=omega0[sx,sy,z0:z1]
        b=omega0[slice(sx.start+dx,sx.stop+dx),slice(sy.start+dy,sy.stop+dy),z0+dz:z1+dz]
        # 2*wbar = weight*(1/omega0i + 1/omega0j).
        accum+=weight*np.sum(1/a+1/b);crossing+=a.size
    E=accum*V*V*Gc/(.25*delta);ratioE=E/(Gc*side**2)
    return dict(n=n,points=n**3,spacing_m=h,horizon_ratio=ratio,horizon_m=delta,side_m=side,Gc_demo_J_per_m2=Gc,crack_area_m2=side**2,crack_energy_J=E,ratio=ratioE,relative_error=ratioE-1,bulk_lattice_ratio=bulk_planar_ratio(ratio),crossing_unordered_bonds=crossing,label='NUMERICAL NORMALIZATION TEST; finite free-edge bias retained')

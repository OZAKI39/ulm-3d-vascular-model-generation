"""Read-only intact BAQP first-variation diagnostics; never a production BC.

Every evaluation requires an explicit material. No phase or equilibrium solver.
"""
from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree
from .kernel import Family
from .energy import constitutive


@dataclass
class Cube:
    family: Family
    ijk: np.ndarray
    L: float
    h: float
    n: int
    face_areas: np.ndarray

    @property
    def distances(self):
        return np.column_stack([self.L/2+self.family.X[:, a]*s
                                for a in range(3) for s in (1, -1)])

    @property
    def full_support(self):
        return self.distances.min(axis=1) >= self.family.horizon


def cube_coordinates(n, L):
    h=L/(n-1);ijk=np.indices((n,n,n)).reshape(3,-1).T
    X=np.linspace(-L/2,L/2,n)[ijk]
    w=np.full(n,h);w[[0,-1]]=h/2
    V=np.prod(w[ijk],axis=1)
    areas=[]
    for a in range(3):
        for side in (0,n-1):
            ids=ijk[:,a]==side
            areas.append(np.prod(w[ijk[ids][:,[b for b in range(3) if b!=a]]],axis=1).sum())
    return X,V,ijk,np.asarray(areas)


def make_cube(n, L, delta):
    X,V,ijk,areas=cube_coordinates(n,L)
    pairs=cKDTree(X).query_pairs(delta,output_type='ndarray')
    return Cube(Family.from_arrays(X,V,pairs,delta),ijk,L,L/(n-1),n,areas)


def continuum(material, L, epsilon):
    F=np.diag([1.,1.,1.+epsilon]);W,P=constitutive(F,material)
    return dict(F=F,E=.5*(F.T@F-np.eye(3)),W_Pa=float(W),P_Pa=P,
                Pzz_Pa=float(P[2,2]),area_m2=L**2,volume_m3=L**3,
                reaction_N=float(P[2,2]*L**2),energy_J=float(W*L**3),
                dU_depsilon_J=float(P[2,2]*L**3))


def nodal(f, pair_force):
    force=np.zeros_like(f.X);i,j=f.pairs.T
    np.add.at(force,i,pair_force);np.add.at(force,j,-pair_force)
    return force


def force_terms(f, state, weight=None):
    """Linear splitting in wbar at fixed intact kinematics and Piola stress.

    Z must be rebuilt for EACH weight component. Splitting only incident direct
    weights would miss the neighbors' contribution through the Z accumulation.
    """
    w=f.wbar if weight is None else np.asarray(weight,float)
    i,j=f.pairs.T;V=f.volume;xi=f.xi;r2=f.length**2
    pair_direct=2*w[:,None]*np.einsum('bij,bj->bi',state.P,xi)/r2[:,None]
    pair_direct*=V[i,None]*V[j,None]
    projection=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r2[:,None,None]
    term=w[:,None,None]*(state.P@projection)
    Z=np.zeros((len(V),3,3))
    np.add.at(Z,i,term*V[j,None,None]);np.add.at(Z,j,term*V[i,None,None])
    pair_z=(np.einsum('bij,bj->bi',Z[i],state.shape_i)
            -np.einsum('bij,bj->bi',Z[j],state.shape_j))*V[i,None]*V[j,None]
    direct=nodal(f,pair_direct);correction=nodal(f,pair_z)
    return dict(direct=direct,correction=correction,total=direct+correction,
                pair_direct=pair_direct,pair_correction=pair_z,Z=Z)


def endpoint_weight(f, mask):
    i,j=f.pairs.T;mask=np.asarray(mask,float)
    return .5*f.omega*(mask[i]/f.omega0[i]+mask[j]/f.omega0[j])


def force_stats(force, mask, R):
    a=np.asarray(force)[mask];mag=np.linalg.norm(a,axis=1)
    if not len(a):return dict(count=0,max_N=None,rms_N=None,max_fz_N=None,max_over_R=None,rms_over_R=None)
    return dict(count=len(a),max_N=float(mag.max()),rms_N=float(np.sqrt(np.mean(mag**2))),
                max_fz_N=float(abs(a[:,2]).max()),max_over_R=float(mag.max()/abs(R)),
                rms_over_R=float(np.sqrt(np.mean(mag**2))/abs(R)))


def required_force(state):
    return -state.force


def self_neighbor_terms(f, state):
    """Split each oriented nodal row by the endpoint supplying its coefficient.

    This is different from splitting an unordered bond by its arbitrary ID order.
    The four nodal arrays sum to the original force without changing any weight.
    """
    i,j=f.pairs.T;V=f.volume
    common=f.omega[:,None]*np.einsum('bij,bj->bi',state.P,f.xi)/f.length[:,None]**2
    common*=V[i,None]*V[j,None]
    di=common/f.omega0[i,None];dj=common/f.omega0[j,None]
    zi=np.einsum('bij,bj->bi',state.Z[i],state.shape_i)*V[i,None]*V[j,None]
    zj=np.einsum('bij,bj->bi',state.Z[j],state.shape_j)*V[i,None]*V[j,None]
    out={k:np.zeros_like(f.X) for k in ('direct_self','direct_neighbor','Z_self','Z_neighbor')}
    for key,a,b in [('direct_self',di,-dj),('direct_neighbor',dj,-di),
                    ('Z_self',zi,zj),('Z_neighbor',-zj,-zi)]:
        np.add.at(out[key],i,a);np.add.at(out[key],j,b)
    return out


def full_lattice_omega0(h, delta):
    """Reference coefficient for algebraic wbar = w_lattice + delta_wbar only."""
    m=int(np.ceil(delta/h));v=np.indices((2*m+1,)*3).reshape(3,-1).T-m
    r=np.linalg.norm(v*h,axis=1);r=r[(r>0)&(r<delta)]
    return float(np.sum((1-r/delta)**2)*h**3)


def cumulative(cube, force, R, side):
    """Signed outward force, integrated over one half (middle plane excluded).

    Quantiles use the observed signed half-reaction, not R_cont. Nonmonotone
    cumulative curves retain their overshoots; quantiles are first crossings,
    never clipped/reweighted into an apparent positive measure.
    """
    z=cube.family.X[:,2];tol=cube.h*1e-8;rows=[];acc=0.
    for layer in sorted(np.unique(cube.ijk[:,2]),reverse=side==1):
        mask=cube.ijk[:,2]==layer;coord=float(z[mask][0])
        if side*coord<=tol:continue
        value=float(side*force[mask,2].sum());acc+=value;depth=cube.L/2-side*coord
        rows.append(dict(side='top' if side==1 else 'bottom',layer=int(layer),depth_m=depth,
                         depth_over_delta=depth/cube.family.horizon,depth_over_h=depth/cube.h,
                         layer_force_N=value,cumulative_N=acc,cumulative_over_R=acc/R))
    quantiles={}
    for q in (.5,.9,.95,.99):
        hits=[r for r in rows if acc>0 and r['cumulative_N']>=q*acc]
        quantiles[f'd{int(100*q)}_m']=hits[0]['depth_m'] if hits else None
    return rows,dict(reaction_N=acc,ratio_to_continuum=acc/R,relative_error=abs(acc/R-1),
                     monotone=all(r['layer_force_N']>=-abs(R)*1e-12 for r in rows),
                     quantile_definition='first layer crossing of fraction of observed signed half-reaction',**quantiles)

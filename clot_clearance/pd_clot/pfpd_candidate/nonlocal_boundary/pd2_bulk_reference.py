"""Selective PD2 oracle: accepted weights, bond law and adjoint functions.

Geometry-only lattice scaffolding follows the old closed-stencil audit. No C1
force, K, Z, bulk coefficients or bulk results are used. Forces away from the
selected observation are partial and must never be used as a full state.
"""
import numpy as np
from ..force_stencil_audit import lattice_offsets,H,DELTA
from ..kernel import radial
from ..higher_order_rk.moment_matrix import local_system
from ..higher_order_rk.rk_kinematics import deformation
from ..higher_order_rk.mechanics_candidate import Reference,GradientOperator
from ..higher_order_rk.variational_force import energy_force
from ..energy import Material
from .nonaffine_manufactured import QuadraticAxial,QuadraticShear
from .local_surrogate import AffineSurrogate


def field_for(kind):
    mat=Material('svk',1000.,2000./3)
    if kind=='nonaffine_axial':return QuadraticAxial()
    if kind in ('nonaffine_shear','classical_shear'):return QuadraticShear()
    F=np.eye(3)
    if kind=='normal_affine':F[2,2]+=.003
    elif kind=='shear_affine':F[0,2]=.003
    else:raise ValueError(kind)
    return AffineSurrogate(F,mat)


class PD2Stencil:
    def __init__(self,n=27,center=(0,0,0)):
        self.n=n;self.ijk=np.indices((n,)*3).reshape(3,-1).T
        self.X=np.linspace(-(n-1)*H/2,(n-1)*H/2,n)[self.ijk]+np.asarray(center)
        self.V=H**3*np.prod(np.where((self.ijk==0)|(self.ijk==n-1),.5,1.),axis=1)
        self.offsets=lattice_offsets();self._neighbors={}

    def id_at_offset(self,offset):return int(np.ravel_multi_index(tuple(np.asarray(offset)+self.n//2),(self.n,)*3))
    def neighbors(self,p):
        if p not in self._neighbors:
            k=self.ijk[p]+self.offsets;k=k[np.all((k>=0)&(k<self.n),axis=1)]
            self._neighbors[p]=np.ravel_multi_index(k.T,(self.n,)*3)
        return self._neighbors[p]

    def dependencies(self,p):
        # Reverse final nodal adjoint: G at p and all sources with target p.
        D0=np.array([p]);D1=np.union1d(D0,self.neighbors(p))
        # Bond adjoint: every incident bond for each required G.
        D2=np.unique(np.concatenate([D1]+[self.neighbors(i) for i in D1]))
        # F at both endpoints and m(endpoint omega0) require complete families.
        D3=np.unique(np.concatenate([D2]+[self.neighbors(i) for i in D2]))
        return D0,D1,D2,D3

    def closure(self,p):
        levels=self.dependencies(p);rows=[]
        for name,ids in zip(['D0','D1','D2','D3'],levels):
            off=self.ijk[ids]-self.ijk[p]
            rows.append(dict(level=name,count=len(ids),maximum_radius_m=float(np.linalg.norm(off*H,axis=1).max()),maximum_lattice_offset=int(abs(off).max())))
        full=all(len(self.neighbors(i))==len(self.offsets) for i in levels[2])
        vol=bool(np.all(self.V[levels[3]]==H**3))
        return dict(N=self.n,observation=self.X[p],levels=rows,D2_families_complete=full,D3_volumes_bulk=vol,
                    gradient_weight_support_complete=full and vol,bond_weight_support_complete=full and vol,
                    all_required_nodes_complete=full and vol,
                    graph='gradU[p] <- direct incident bonds and G[D1] <- P[bonds incident D1], g[D1] <- F[D2], m[D2] <- positions/volumes[D3]')

    def force(self,p,field):
        D0,D1,D2,D3=self.dependencies(p);X=self.X;V=self.V
        sources=[];targets=[];weights=[];omega0={}
        for i in D2:
            j=self.neighbors(i);xi=X[j]-X[i];w=radial(np.linalg.norm(xi,axis=1)/DELTA)
            _,_,rank,_,_,g=local_system(xi,w,V[j],DELTA,9)
            if rank!=9:raise ValueError('PD2 local rank lost')
            sources.append(np.full(len(j),i));targets.append(j);weights.append(g)
            # Same family convention; summation-rounding is checked against full evaluator.
            omega0[i]=float(np.sum(w*V[j]))
        source=np.concatenate(sources);target=np.concatenate(targets);g=np.concatenate(weights)
        x=field.position(X);F=deformation(X,x,source,target,g)
        pairs=np.unique(np.concatenate([np.sort(np.column_stack([np.full(len(self.neighbors(i)),i),self.neighbors(i)]),axis=1) for i in D1]),axis=0)
        i,j=pairs.T;xi=X[j]-X[i];r2=np.sum(xi*xi,axis=1)
        m=radial(np.sqrt(r2)/DELTA)*(np.array([1/omega0[k] for k in i])+np.array([1/omega0[k] for k in j]))*V[i]*V[j]
        ref=Reference(X,V,pairs,xi,r2,m)
        needed=np.isin(source,D1);op=GradientOperator('PD2',source[needed],target[needed],g[needed])
        _,force=energy_force(ref,op,x,F,Material('svk',1000.,2000./3))
        return force[p]


def oracle_job(spec):
    offset=np.asarray(spec.get('offset',[0,0,0]));point=np.asarray(spec.get('point',[0,0,spec.get('Z',0.)]))
    s=PD2Stencil(spec.get('n',27),point-offset*H);p=s.id_at_offset(offset)
    force=s.force(p,field_for(spec['kind']))
    return dict(**spec,force_N=force.tolist(),q_bulk_N_m3=(force/s.V[p]).tolist(),closure=s.closure(p))

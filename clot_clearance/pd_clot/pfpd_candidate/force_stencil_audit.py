"""Selective, intact affine BAQP diagnostics; not a production or BC evaluator.

Input coordinates are affine. Local Fbar/Ftilde/P reconstruction is retained,
including floating-point effects, to match the original evaluate() path.
The kernel, endpoint volumes, unordered bonds, K and Z are unchanged.
No equilibrium, phase, fracture, or global N19/N27 bond state is constructed.
"""
from functools import lru_cache
import numpy as np
from .energy import Material, constitutive
from .kernel import radial

H = 0.0000625
DELTA = 0.00025125
EPSILON = 0.003
FORCE_TOL = 1e-10
COEFFICIENT_TOL = 1e-12
RECONSTRUCTION_TOL = 1e-12


def material():
    return Material(model='svk', mu=1000., bulk=666.6666666666666)


def configuration():
    return dict(task='PFPD_BAQP_FULLY_CLOSED_FORCE_STENCIL_AUDIT_V1',
                model='PFPD_AT1_SPECTRAL_CANDIDATE_V1', physics='INTACT_MECHANICS_ONLY', phase='DISABLED', s=0,
                h_m=H, delta_m=DELTA, delta_over_h=DELTA/H, engineering_strain=EPSILON,
                mu_Pa=1000., lambda_Pa=0., bulk_Pa=material().bulk, kernel='quadratic',
                grids=[11,19,27], local_normalization='abs(Pzz) * h**2',
                primary_force_tolerance=FORCE_TOL, coefficient_match_tolerance=COEFFICIENT_TOL,
                decomposition_reconstruction_tolerance=RECONSTRUCTION_TOL,
                full_pytest='DO_NOT_RUN', GPU='NOT_USED_BY_DESIGN', CPU_workers=1, BLAS_threads=1,
                off_center_selection='Nearest lattice depths round(k*delta/h), k=1,2,3; plus ceil(delta/h), ceil(2*delta/h) to test asymmetric closed one/two-hop stencils; reflected pairs.',
                center_cancellation_caveat='A zero center force does not establish translation-invariant coefficients or minimum sufficient graph closure.')


@lru_cache(maxsize=8)
def lattice_offsets(h=H, delta=DELTA):
    m=int(np.ceil(delta/h))
    a=np.indices((2*m+1,)*3).reshape(3,-1).T-m
    radius=np.linalg.norm(a*h,axis=1)
    # Same <=delta selection as cKDTree.query_pairs; radial() is zero at delta.
    a=a[(radius>0)&(radius<=delta)]
    a.setflags(write=False)
    return a


@lru_cache(maxsize=8)
def dependency_offsets(h=H, delta=DELTA):
    """Exact-length walks Dq=N(Dq-1), not disjoint distance shells."""
    offsets=lattice_offsets(h,delta)
    levels=[np.zeros((1,3),dtype=int)]
    for _ in range(3):
        levels.append(np.unique((levels[-1][:,None,:]+offsets[None,:,:]).reshape(-1,3),axis=0))
    for a in levels:a.setflags(write=False)
    return tuple(levels)


def bulk_oracle(epsilon=EPSILON):
    offsets=lattice_offsets();xi=offsets*H;r=np.linalg.norm(xi,axis=1);w=radial(r/DELTA)
    V=H**3;omega0=float(np.sum(w*V))
    K=np.einsum('n,ni,nj->ij',w*V,xi,xi);inv=np.linalg.inv(K)
    F=np.diag([1.,1.,1.+epsilon]);_,P=constitutive(F,material())
    wb=w/omega0;proj=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r[:,None,None]**2
    Z=np.sum(wb[:,None,None]*(P@proj)*V,axis=0)
    B=w[:,None]*(xi@inv.T)
    # Each incident unordered physical bond contributes once to the center.
    direct=2*wb[:,None]*(xi@P.T)/r[:,None]**2*V**2
    correction=2*(B@Z.T)*V**2
    force=np.sum(direct+correction,axis=0)
    _,reference_P=constitutive(np.diag([1.,1.,1.+EPSILON]),material())
    local_scale=abs(reference_P[2,2])*H**2
    return dict(neighbor_offset_count=len(offsets),offsets=offsets,xi=xi,omega=w,wbar=wb,
                omega0_inf=omega0,K_inf=K,K_inverse_inf=inv,K_eigenvalues=np.linalg.eigvalsh(K),
                Z_inf=Z,F=F,P=P,Pzz=float(P[2,2]),local_force_scale_N=float(local_scale),
                center_force_N=force,center_force_norm_N=float(np.linalg.norm(force)),
                center_force_local_normalized=float(np.linalg.norm(force)/local_scale),
                direct_N=direct.sum(0),Z_correction_N=correction.sum(0),
                offset_inclusion='0 < norm(offset*h) <= delta; radial quadratic has rho<1 support',
                PASS=bool(np.linalg.norm(force)/local_scale<=FORCE_TOL))


class AffineStencil:
    def __init__(self,n,epsilon=EPSILON):
        if n<3 or n%2!=1:raise ValueError('Odd boundary-node cube required')
        self.n=n;self.h=H;self.delta=DELTA;self.L=(n-1)*H;self.epsilon=epsilon
        self.ijk=np.indices((n,n,n)).reshape(3,-1).T
        self.X=np.linspace(-self.L/2,self.L/2,n)[self.ijk]
        weights=np.full(n,H);weights[[0,-1]]=H/2
        self.V=np.prod(weights[self.ijk],axis=1)
        if abs(self.V.sum()/self.L**3-1)>2e-15:raise ValueError('Volume quadrature mismatch')
        self.offsets=lattice_offsets();self.oracle=bulk_oracle(epsilon)
        self.P=self.oracle['P'];self.F=self.oracle['F'];self.scale=self.oracle['local_force_scale_N']
        self.x=self.X@self.F.T;self.x[:,2]+=epsilon*self.L/2
        self.u=self.x-self.X
        self._neighbors={};self._omega={};self._K={};self._inv={};self._Z={};self._Fbar={}

    def id_at_offset(self,offset):
        index=np.asarray(offset,dtype=int)+self.n//2
        if np.any(index<0) or np.any(index>=self.n):raise ValueError('Point outside cube')
        return int(np.ravel_multi_index(tuple(index),(self.n,)*3))

    def ids(self,ijk):
        return np.ravel_multi_index(np.asarray(ijk).T,(self.n,)*3)

    def neighbors(self,pid):
        if pid not in self._neighbors:
            candidate=self.ijk[pid]+self.offsets
            available=np.all((candidate>=0)&(candidate<self.n),axis=1)
            ids=self.ids(candidate[available])
            # Check real floating-coordinate inclusion as in Family construction.
            ids=ids[np.linalg.norm(self.X[ids]-self.X[pid],axis=1)<=self.delta]
            self._neighbors[pid]=ids
        return self._neighbors[pid]

    def row(self,pid):
        ids=self.neighbors(pid)
        # Lexicographic i<j bonds: first the lower-endpoint assembly pass,
        # then the upper-endpoint pass. No global bond list is needed.
        ids=np.r_[ids[ids>pid],ids[ids<pid]]
        xi=self.X[ids]-self.X[pid]
        radius=np.linalg.norm(xi,axis=1)
        return ids,xi,radius,radial(radius/self.delta)

    def omega0(self,pid):
        if pid not in self._omega:
            ids,_,_,w=self.row(pid);self._omega[pid]=float(np.add.accumulate(w*self.V[ids])[-1])
        return self._omega[pid]

    def K(self,pid):
        if pid not in self._K:
            ids,xi,_,w=self.row(pid)
            terms=w[:,None,None]*xi[:,:,None]*xi[:,None,:]*self.V[ids,None,None]
            K=np.add.accumulate(terms,axis=0)[-1]
            eig=np.linalg.eigvalsh(K)
            if eig[0]<=max(eig[-1],0)*1e-12:raise ValueError('Rank deficient family')
            self._K[pid]=K;self._inv[pid]=np.linalg.inv(K)
        return self._K[pid]

    def inverse(self,pid):
        self.K(pid)
        return self._inv[pid]

    def Fbar(self,pid):
        if pid not in self._Fbar:
            ids,xi,_,w=self.row(pid)
            B=w[:,None]*np.einsum('ij,nj->ni',self.inverse(pid),xi)
            du=self.u[ids]-self.u[pid]
            terms=du[:,:,None]*B[:,None,:]*self.V[ids,None,None]
            F=np.eye(3)[None,:,:].copy()
            np.add.at(F,np.zeros(len(ids),dtype=int),terms)
            self._Fbar[pid]=F[0]
        return self._Fbar[pid]

    def bond_stress(self,low,high):
        """Small incident-bond batch; never an all-domain P/Ftilde array."""
        Fi=np.array([self.Fbar(int(j)) for j in low]);Fj=np.array([self.Fbar(int(j)) for j in high])
        av=.5*(Fi+Fj);xi=self.X[high]-self.X[low];eta=self.x[high]-self.x[low]
        r2=np.linalg.norm(xi,axis=1)**2
        residual=eta-np.einsum('bij,bj->bi',av,xi)
        Ft=av+residual[:,:,None]*xi[:,None,:]/r2[:,None,None]
        return constitutive(Ft,material())[1]

    def weight(self,pid,ids,w,channel):
        if channel=='bulk_weight':return w/self.oracle['omega0_inf']
        actual=.5*w*(1/self.omega0(pid)+1/np.array([self.omega0(int(j)) for j in ids]))
        if channel=='original':return actual
        if channel=='normalization_excess':return actual-w/self.oracle['omega0_inf']
        raise ValueError(channel)

    def Z(self,pid,channel='original'):
        key=(pid,channel)
        if key not in self._Z:
            ids,xi,r,w=self.row(pid)
            projection=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r[:,None,None]**2
            P=self.bond_stress(np.minimum(pid,ids),np.maximum(pid,ids))
            for component in ('original','bulk_weight','normalization_excess'):
                wb=self.weight(pid,ids,w,component)
                terms=wb[:,None,None]*(P@projection)*self.V[ids,None,None]
                self._Z[(pid,component)]=np.add.accumulate(terms,axis=0)[-1]
        return self._Z[key]

    def coefficient_sets(self,pid):
        first=np.union1d([pid],self.neighbors(pid)).astype(int)
        second=np.unique(np.concatenate([first]+[self.neighbors(int(j)) for j in first]))
        return first,second

    def _force_channel(self,pid,channel):
        neighbors,_,r,w=self.row(pid)
        # Same unordered physical-bond convention and ti-tj assembly as evaluate.
        low=np.minimum(pid,neighbors);high=np.maximum(pid,neighbors)
        xi=self.X[high]-self.X[low];sign=np.where(low==pid,1.,-1.)
        wb=self.weight(pid,neighbors,w,channel)
        inv_low=np.array([self.inverse(int(j)) for j in low]);inv_high=np.array([self.inverse(int(j)) for j in high])
        Bi=w[:,None]*np.einsum('nij,nj->ni',inv_low,xi)
        Bj=-w[:,None]*np.einsum('nij,nj->ni',inv_high,xi)
        Zi=np.array([self.Z(int(j),channel) for j in low]);Zj=np.array([self.Z(int(j),channel) for j in high])
        P=self.bond_stress(low,high)
        direct=wb[:,None]*np.einsum('bij,bj->bi',P,xi)/r[:,None]**2
        zi=np.einsum('nij,nj->ni',Zi,Bi);zj=np.einsum('nij,nj->ni',Zj,Bj)
        signed_volume=(sign*self.V[low]*self.V[high])[:,None]
        force=np.sum(((direct+zi)-(-direct+zj))*signed_volume,axis=0)
        return dict(force=force,direct=np.sum(2*direct*signed_volume,axis=0),
                    Z_correction=np.sum((zi-zj)*signed_volume,axis=0))

    def evaluate_point(self,pid):
        channels={k:self._force_channel(pid,k) for k in ('original','bulk_weight','normalization_excess')}
        primary=channels['original'];f=primary['force']
        errors={k:float(np.linalg.norm(v['force']-v['direct']-v['Z_correction'])/self.scale) for k,v in channels.items()}
        errors['weight_split']=float(np.linalg.norm(f-channels['bulk_weight']['force']-channels['normalization_excess']['force'])/self.scale)
        if max(errors.values())>RECONSTRUCTION_TOL:raise ValueError('Decomposition reconstruction failed')
        depth=int(np.minimum(self.ijk[pid],self.n-1-self.ijk[pid]).min())
        return dict(ID=int(pid),ijk=self.ijk[pid],X_m=self.X[pid],boundary_distance_cells=depth,
                    boundary_distance_m=depth*H,boundary_distance_over_delta=depth*H/DELTA,
                    force_N=f,force_norm_N=float(np.linalg.norm(f)),
                    force_local_normalized=float(np.linalg.norm(f)/self.scale),
                    abs_fz_local_normalized=float(abs(f[2])/self.scale),
                    channels=channels,decomposition_errors_local_scaled=errors)

    def closure(self,pid):
        theoretical=dependency_offsets();origin=self.ijk[pid]
        actual=np.array([origin]);rows=[]
        for depth,relative in enumerate(theoretical):
            required=origin+relative
            exists=np.all((required>=0)&(required<self.n),axis=1)
            ids=self.ids(required[exists]);inside_volume=self.V[ids]
            if depth:
                candidate=(actual[:,None,:]+self.offsets[None,:,:]).reshape(-1,3)
                candidate=candidate[np.all((candidate>=0)&(candidate<self.n),axis=1)]
                actual=np.unique(candidate,axis=0)
            all_exist=bool(exists.all())
            all_interior=bool(all_exist and np.all(inside_volume==H**3))
            radius=np.linalg.norm(relative*H,axis=1)
            rows.append(dict(depth=depth,required_unique_nodes=len(relative),finite_graph_unique_nodes=len(actual),
                             missing_required_nodes=int((~exists).sum()),max_required_radius_m=float(radius.max()),
                             max_required_coordinate_offset_cells=int(abs(relative).max()),
                             max_finite_graph_radius_m=float(np.linalg.norm((actual-origin)*H,axis=1).max()),
                             all_required_nodes_exist=all_exist,all_required_volumes_are_interior=all_interior,
                             nonbulk_volume_node_count=int(np.count_nonzero(inside_volume!=H**3)),
                             closed_with_bulk_volumes=all_interior))
        required_D2=origin+np.unique(np.vstack(theoretical[:3]),axis=0)
        valid=np.all((required_D2>=0)&(required_D2<self.n),axis=1)
        family_failures=[]
        for node in required_D2[valid]:
            nid=int(self.ids(node[None,:])[0])
            if not np.array_equal(self.ijk[self.neighbors(nid)]-node,self.offsets):family_failures.append(nid)
        D2_full=bool(valid.all() and not family_failures)
        fully=bool(rows[3]['closed_with_bulk_volumes'] and D2_full)
        return dict(ID=int(pid),center_ijk=origin,levels=rows,
                    definition='D0={i}; Dq=union of family neighbors of D(q-1). Sets overlap; counts are not disjoint shells. Required sets use the infinite lattice.',
                    all_required_D2_nodes_present=bool(valid.all()),every_D2_node_has_complete_family=D2_full,
                    present_D2_incomplete_family_count=len(family_failures),actual_required_graph_depth=3,
                    required_volumes='Every required D0..D3 node must have h^3, not boundary quadrature weights.',
                    FULL_FORCE_STENCIL_CLOSED=fully)

    def coefficients(self,pid):
        first,second=self.coefficient_sets(pid);o=self.oracle
        omega=np.array([self.omega0(int(j)) for j in second]);K=np.array([self.K(int(j)) for j in first]);Z=np.array([self.Z(int(j)) for j in first])
        errors=dict(omega0=float(abs(omega/o['omega0_inf']-1).max()),
                    K=float(np.linalg.norm(K-o['K_inf'],axis=(1,2)).max()/np.linalg.norm(o['K_inf'])),
                    Z=float(np.linalg.norm(Z-o['Z_inf'],axis=(1,2)).max()/max(np.linalg.norm(o['Z_inf']),np.finfo(float).tiny)))
        wbar_error=0.
        for j in first:
            ids,_,_,w=self.row(int(j));actual=self.weight(int(j),ids,w,'original');ref=w/o['omega0_inf']
            wbar_error=max(wbar_error,float(np.max(abs(actual/ref-1))))
        errors['wbar']=wbar_error
        # Center and the six extremal D1 offsets are chosen by geometry, not IDs.
        representatives=[pid]
        for axis in range(3):
            candidates=self.neighbors(pid)
            for choose in (np.argmin,np.argmax):representatives.append(int(candidates[choose(self.X[candidates,axis])]))
        details=[]
        for j in dict.fromkeys(representatives):
            Kj=self.K(j);Zj=self.Z(j)
            details.append(dict(ID=int(j),X_m=self.X[j],omega0_ratio=self.omega0(j)/o['omega0_inf'],
                                K=Kj,Z=Zj,K_relative_error=float(np.linalg.norm(Kj-o['K_inf'])/np.linalg.norm(o['K_inf'])),
                                Z_relative_error=float(np.linalg.norm(Zj-o['Z_inf'])/max(np.linalg.norm(o['Z_inf']),np.finfo(float).tiny))))
        K_kinematic_error=float(max(np.linalg.norm(self.K(int(j))-o['K_inf'])/np.linalg.norm(o['K_inf']) for j in second))
        errors['K_direct_B']=errors['K']
        errors['K']=max(errors['K'],K_kinematic_error)
        return dict(max_relative_error=errors,omega0_required_nodes=len(second),K_Z_required_nodes=len(first),
                    kinematic_K_D2_max_relative_error=K_kinematic_error,
                    omega0_D0_D1_max_error=float(max(abs(self.omega0(int(j))/o['omega0_inf']-1) for j in first)),
                    scope='omega0 and K: D0 union D1 union D2 (K includes Fbar reconstruction); Z: D0 union D1; K_direct_B separately restricts K to D0 union D1; wbar: all bonds incident to D0 union D1. Matrix errors are relative Frobenius norms.',
                    representative_coefficients=details)

    def stress_roundoff_diagnostic(self,pid):
        """Supplementary constant-P coefficients; never replaces the strict gate.

        This isolates the difference introduced by reconstructing a prescribed
        affine F in Float64. Force and primary coefficient gates retain the
        original reconstructed-P values.
        """
        first,second=self.coefficient_sets(pid);ref=self.oracle;rows=[]
        for j in first:
            ids,xi,r,w=self.row(int(j));wb=self.weight(int(j),ids,w,'original')
            proj=np.eye(3)-xi[:,:,None]*xi[:,None,:]/r[:,None,None]**2
            Z_constant=np.add.accumulate(wb[:,None,None]*(self.P@proj)*self.V[ids,None,None],axis=0)[-1]
            rows.append((np.linalg.norm(Z_constant-ref['Z_inf'])/np.linalg.norm(ref['Z_inf']),
                         np.linalg.norm(self.Z(int(j))-Z_constant)/np.linalg.norm(ref['Z_inf'])))
        return dict(constant_P_Z_max_bulk_relative_error=float(max(r[0] for r in rows)),
                    reconstructed_P_minus_constant_P_Z_max_relative_error=float(max(r[1] for r in rows)),
                    Fbar_D2_max_absolute_error=float(max(np.max(abs(self.Fbar(int(j))-self.F)) for j in second)),
                    primary_gate_replaced=False,interpretation='The constant-P result is a roundoff diagnostic only; the original reconstructed-P coefficient gate is retained.')

    def cache_counts(self):
        return dict(neighbor_rows=len(self._neighbors),omega0_nodes=len(self._omega),K_nodes=len(self._K),
                    Fbar_nodes=len(self._Fbar),Z_tensors=len(self._Z),full_domain_bond_state_allocated=False)


def validate_selective(n):
    """Small full-state gate, including all coefficients needed by chosen points."""
    from scipy.spatial import cKDTree
    from .kernel import Family
    from .baqp import evaluate
    from .boundary_force_diagnostics import force_terms
    grid=AffineStencil(n)
    pairs=cKDTree(grid.X).query_pairs(DELTA,output_type='ndarray')
    pairs=pairs[np.lexsort((pairs[:,1],pairs[:,0]))]
    full=Family.from_arrays(grid.X,grid.V,pairs,DELTA)
    x=grid.X@grid.F.T;x[:,2]+=EPSILON*grid.L/2
    state=evaluate(full,x,material());terms=force_terms(full,state)
    offsets=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,-1,-1)]
    rows=[];first_all=set();second_all=set()
    for offset in offsets:
        pid=grid.id_at_offset(offset);value=grid.evaluate_point(pid);f=value['channels']['original']
        errors={name:float(np.linalg.norm(f[name]-reference)/grid.scale) for name,reference in
                [('force',state.force[pid]),('direct',terms['direct'][pid]),('Z_correction',terms['correction'][pid])]}
        first,second=grid.coefficient_sets(pid);first_all.update(map(int,first));second_all.update(map(int,second))
        rows.append(dict(offset_cells=offset,ID=pid,full_force_N=state.force[pid],selective_force_N=f['force'],local_scaled_errors=errors))
    o=grid.oracle
    coefficients=dict(omega0=max(abs(grid.omega0(j)-full.omega0[j])/o['omega0_inf'] for j in second_all),
                      K=max(np.linalg.norm(grid.K(j)-state.K[j])/np.linalg.norm(o['K_inf']) for j in first_all),
                      Z=max(np.linalg.norm(grid.Z(j)-state.Z[j])/np.linalg.norm(o['Z_inf']) for j in first_all))
    force_max=max(max(r['local_scaled_errors'].values()) for r in rows)
    passed=force_max<=RECONSTRUCTION_TOL and max(coefficients.values())<=COEFFICIENT_TOL
    return dict(n=n,particles=n**3,full_bonds=len(pairs),observation_points=rows,
                force_max_local_scaled_error=force_max,coefficient_max_scaled_errors=coefficients,
                coefficient_nodes_checked=dict(omega0=len(second_all),K_and_Z=len(first_all)),
                force_error_scale_N=grid.scale,coefficient_error_scales='Infinite-lattice omega0 and Frobenius norms of K_inf, Z_inf',
                affine_specialization='Affine input coordinates; local Fbar/Ftilde/P reconstructed with the existing formulas. Both evaluators use deterministic lexicographic unordered bonds.',
                PASS=bool(passed))

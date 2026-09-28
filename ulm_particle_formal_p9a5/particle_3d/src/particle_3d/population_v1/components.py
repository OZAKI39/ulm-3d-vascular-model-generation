"""Canonical connected-component resistance solve; no global dense 6N array."""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu, eigsh, LinearOperator, lsmr
from .. import particle65_motion as p65
from ..particle6_stepper import bind_query_dependency
from ..particle9a_motion import augment_planar_system
from ..resistance_assembly import ResistanceSystem
from ..resistance_solver import solve_resistance, ResistanceSolution, ResistanceSystemIllConditioned
from ..hydrodynamic_resistance import sphere_self_diagonal
from ..nearfield_handoff import exact_interactions
from ..wall_gap import wall_gap
from .pairs import candidates, components, sphere_gap


@dataclass
class SparseRow:
    indices: np.ndarray
    values: np.ndarray
    def __matmul__(self, vector):return float(self.values@vector[self.indices])


def sparse_assembly(shapes, free, mu, wall, *, policy):
    ids=tuple(sorted(shapes));index={i:6*k for k,i in enumerate(ids)}
    interactions=bind_query_dependency(exact_interactions,dict(all_pairs=candidates,pair_gap=sphere_gap))
    specs,_,_=interactions(shapes,wall,policy,mu)
    diag=np.concatenate([sphere_self_diagonal(shapes[i],mu) for i in ids])
    u=np.concatenate([free[i] for i in ids]);rr=[];cc=[];vv=[];blocks=[];rows=[]
    for spec in specs:
        z,n,record=policy.evaluate(spec,shapes,mu)
        ix=np.arange(index[spec.particle_i_id],index[spec.particle_i_id]+3);v=n
        if spec.particle_j_id is not None:
            ix=np.r_[ix,np.arange(index[spec.particle_j_id],index[spec.particle_j_id]+3)];v=np.r_[n,-n]
        if z:
            rr.extend(np.repeat(ix,len(ix)));cc.extend(np.tile(ix,len(ix)));vv.extend((z*np.outer(v,v)).ravel())
        blocks.append(record);rows.append(SparseRow(ix,v))
    matrix=sparse.diags(diag,format='csr')+sparse.coo_matrix((vv,(rr,cc)),shape=(len(diag),len(diag))).tocsr()
    return ResistanceSystem(ids,matrix,diag*u,diag,u,blocks,rows)


def sparse_jacobian(system, particles, constraints):
    ordered=sorted(constraints,key=lambda c:c.canonical_id);index={i:6*k for k,i in enumerate(system.ids)}
    rr=[];cc=[];vv=[]
    for k,c in enumerate(ordered):
        if c.offset_m_s!=0:raise ValueError('Only stationary normal constraints supported')
        for pid,sign in [(c.particle_i_id,1),(c.particle_j_id,-1)]:
            if pid is None:continue
            rr.extend([k]*3);cc.extend(range(index[pid],index[pid]+3));vv.extend(sign*c.normal)
    return sparse.coo_matrix((vv,(rr,cc)),shape=(len(ordered),len(system.free))).tocsr(),ordered


def solve_large(system, *, particles=None, constraints=()):
    """Sparse LU plus matrix-free contact Schur active set.

    No dense R, J, Gram, eigenvector basis or inverse is constructed. Redundant
    active rows are retained and solved by LSMR's minimum-norm least squares;
    all original rows are checked at the end, including opposing normals.
    Large-cluster tolerances are explicit FP64 backward-error budgets.
    """
    eps=np.finfo(float).eps;n=len(system.free);d=1/np.sqrt(system.self_diagonal)
    a=(sparse.diags(d)@system.matrix@sparse.diags(d)).tocsc();rhs=d*system.rhs
    sym=float(np.max(np.abs((a-a.T).data),initial=0.))
    norm=float(np.max(np.asarray(abs(a).sum(axis=1))))
    if sym>64*eps*max(1,norm):raise ValueError('NONSYMMETRIC_COMPONENT')
    lo=float(eigsh(a,k=1,which='SA',return_eigenvectors=False,tol=1e-12,v0=np.ones(n))[0])
    hi=float(eigsh(a,k=1,which='LA',return_eigenvectors=False,tol=1e-12,v0=np.ones(n))[0])
    condition=hi/lo if lo>0 else float('inf')
    if not np.isfinite(condition) or eps*n*condition>=1/512:
        raise ResistanceSystemIllConditioned(dict(condition_estimate=condition,solver='SPARSE_COMPONENT'))
    lu=splu(a,permc_spec='NATURAL');hydro=d*lu.solve(rhs);u=hydro.copy()
    scale=max(np.max(np.abs(system.free.reshape(-1,6)[:,:3])),np.max(np.abs(hydro.reshape(-1,6)[:,:3])),np.finfo(float).tiny)
    j,ordered=sparse_jacobian(system,particles,constraints);m=len(ordered)
    budget=512*eps*max(n,m,1)*max(condition,1.)*scale
    lam=np.zeros(m);active=[];iterations=0
    def response(v):return d*lu.solve(d*np.asarray(j.T@v).ravel())
    if m:
        h=np.asarray(j@hydro).ravel()
        diag_bound=max(float(np.max(np.asarray(j.multiply(j)@ (d*d)))),np.finfo(float).tiny)/lo
        dual_budget=budget/diag_bound
        for iterations in range(16*(m+1)**2):
            speeds=np.asarray(j@u).ravel()
            bad=[k for k in range(m) if k not in active and speeds[k]<-budget]
            if not bad:break
            active.append(min(bad,key=lambda k:(speeds[k],k)));active.sort()
            for _ in range(16*(m+1)**2):
                ix=np.asarray(active,dtype=int)
                def product(q):
                    full=np.zeros(m);full[ix]=q
                    return np.asarray(j@response(full)).ravel()[ix]
                op=LinearOperator((len(ix),len(ix)),matvec=product,rmatvec=product,dtype=float)
                candidate=np.zeros(m)
                candidate[ix]=lsmr(op,-h[ix],atol=8*eps,btol=8*eps,conlim=1/eps,maxiter=max(100,8*len(ix)))[0]
                if np.max(np.abs(product(candidate[ix])+h[ix]),initial=0)>budget:
                    raise ResistanceSystemIllConditioned(dict(reason='SPARSE_ACTIVE_KKT_RESIDUAL',active=active))
                negative=[k for k in active if candidate[k]<-dual_budget]
                if not negative:lam=np.maximum(candidate,0.);break
                alpha=min(lam[k]/(lam[k]-candidate[k]) for k in negative)
                lam+=alpha*(candidate-lam)
                leave=[k for k in active if lam[k]<=dual_budget and candidate[k]<0]
                if not leave:raise ResistanceSystemIllConditioned(dict(reason='SPARSE_ACTIVE_SET_STALLED'))
                for k in leave:active.remove(k);lam[k]=0.
            else:raise ResistanceSystemIllConditioned(dict(reason='SPARSE_INNER_ITERATION_LIMIT'))
            u=hydro+response(lam)
        else:raise ResistanceSystemIllConditioned(dict(reason='SPARSE_OUTER_ITERATION_LIMIT'))
    speeds=np.asarray(j@u).ravel();stationarity=a@(u/d)-rhs-d*np.asarray(j.T@lam).ravel()
    backward=float(np.linalg.norm(stationarity)/(norm*np.linalg.norm(u/d)+np.linalg.norm(rhs)+np.finfo(float).tiny))
    comp=float(np.max(np.abs(lam*speeds),initial=0.))
    if np.min(speeds,initial=0)<-budget or comp>budget*max(scale,np.max(np.abs(lam),initial=0.)) or backward>512*eps*n:
        raise ResistanceSystemIllConditioned(dict(reason='SPARSE_FINAL_KKT_FAILED',backward=backward,complementarity=comp))
    record=dict(solver='SPARSE_LU_MATRIX_FREE_SCHUR_ACTIVE_SET',condition_estimate=condition,
        relative_backward_residual=backward,contact_count=m,contact_ids=[c.canonical_id for c in ordered],
        multipliers=lam.tolist(),normal_speeds_m_s=speeds.tolist(),
        contact_kkt=dict(velocity_budget_m_s=budget,iterations=iterations+1),
        max_constraint_violation_m_s=max(0.,-float(np.min(speeds,initial=0.))),
        max_complementarity_residual=comp,global_dense_conversion=False,dissipation=system.dissipation(u))
    return ResistanceSolution(u,hydro,record)


@dataclass
class ComponentSystem:
    ids: tuple
    free: np.ndarray
    blocks: list
    systems: list


def assemble_components(shapes,free,mu,wall,*,policy,gradient_at,raw_candidates=None):
    ids=tuple(sorted(shapes));systems=[];blocks=[]
    pairs=candidates(shapes)
    # Extra broadphase edges are safe, and include contact shells for very small radii.
    for group in components(ids,pairs):
        local={i:shapes[i] for i in group};lf={i:free[i] for i in group}
        if len(group)<=32:
            interactions=bind_query_dependency(exact_interactions,dict(all_pairs=candidates,pair_gap=sphere_gap))
            build=bind_query_dependency(p65.assemble_v1,dict(exact_interactions=interactions))
        else:build=sparse_assembly
        system=build(local,lf,mu,wall,policy=policy)
        system=augment_planar_system(system,local,mu,wall,gradient_at,wall_gap)
        systems.append(system);blocks.extend(system.blocks)
    return ComponentSystem(ids,np.concatenate([free[i] for i in ids]),blocks,systems)


def solve_components(system, *, particles=None, constraints=()):
    index={i:k for k,i in enumerate(system.ids)};u=np.empty_like(system.free);hydro=np.empty_like(u);audits=[]
    for local in system.systems:
        ids=set(local.ids);cs=[c for c in constraints if c.particle_i_id in ids]
        if any(c.particle_j_id is not None and c.particle_j_id not in ids for c in cs):
            raise ValueError('CONTACT_EDGE_MISSING_FROM_COMPONENT_GRAPH')
        solved=(solve_resistance if len(ids)<=32 else solve_large)(local,particles=particles,constraints=cs)
        for k,i in enumerate(local.ids):
            u[6*index[i]:6*index[i]+6]=solved.velocity[6*k:6*k+6]
            hydro[6*index[i]:6*index[i]+6]=solved.unconstrained[6*k:6*k+6]
        audits.append(dict(particle_ids=local.ids,**solved.record))
    return ResistanceSolution(u,hydro,dict(solver='CANONICAL_INTERACTION_COMPONENTS',components=audits,
        component_count=len(audits),maximum_component_size=max(map(lambda x:len(x.ids),system.systems)),
        global_dense_conversion=False,contact_count=sum(a['contact_count'] for a in audits)))

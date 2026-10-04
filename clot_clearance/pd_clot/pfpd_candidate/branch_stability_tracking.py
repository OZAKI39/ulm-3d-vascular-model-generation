"""Stationary continuation diagnostics built on the unchanged coupling audit.

No energy, force, projected KKT or Hessian implementation is duplicated here.
The only new derivative is the chain rule for the stationary residual.
"""
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import eigh,solve
from .coupling_stability_audit import (ROOT,AuditProblem,read_json,save_json,read_npz,
    physics_config,audit_config,projected_kkt,green_minimum)
from scripts import run_monolithic_stability_audit as legacy_hessian

OUT=ROOT/'verification/pfpd_at1_branch_stability_tracking_v1'


def settings():return read_json(OUT/'provenance/PREDECLARED_CONFIG.json')


class BranchDomainExit(RuntimeError):
    def __init__(self,reason,z,metrics):
        super().__init__(reason);self.z=np.array(z,copy=True);self.metrics=metrics


def residual_components(p,z):
    """Same density-normalized projected KKT as the previous audit."""
    x,s=p.unpack(z);state=p.state(x,s);_,g=p.objective_gradient(z)
    q=g[p.nx:]*p.energy_scale/p.weights
    qs=q/np.maximum(state.W0,p.yc)
    return state.force[p.free].ravel()/p.config['mechanical_equilibrium']['force_floor_N'],projected_kkt(s,qs,p.lower,p.sc)


def stationary_residual(p,z,check_domain=True):
    if check_domain:
        x,s=p.unpack(z);m=p.metrics(x,s)
        if m['maximum_s']>=p.sc:raise BranchDomainExit('HIGH_DAMAGE_GATE_REACHED',z,m)
        if m['minimum_Green_strain']<p.config['tolerances']['tensile_min_green_strain']:
            raise BranchDomainExit('PURE_TENSION_DOMAIN_VIOLATION',z,m)
        if m['minimum_lower_margin']<-settings()['bounds_tolerance']:
            raise BranchDomainExit('IRREVERSIBILITY_BOUND_VIOLATION',z,m)
    return np.concatenate(residual_components(p,z))


def residual_jacobian(p,z):
    """Chain rule only; raw energy Hessian comes from the existing helper.

    Algebraic derivative probes are not accepted trial states. At the projected
    residual kink the interior derivative is selected, without claiming ordinary
    differentiability there. SciPy solves residual equations, not energy minima.
    """
    H,_=legacy_hessian.full_hessian(p,z,settings()['residual_jacobian_epsilon'])
    x,s=p.unpack(z);state=p.state(x,s);_,g=p.objective_gradient(z)
    q=g[p.nx:]*p.energy_scale/p.weights;Y=state.W0;D=np.maximum(Y,p.yc)
    dq=H[p.nx:]*p.energy_scale/p.weights[:,None]
    dY=-dq/(2*(1-s[:,None]));dY[:,p.nx:]+=np.diag(Y/(1-s))
    dD=dY*(Y>p.yc)[:,None]
    dqs=dq/D[:,None]-q[:,None]*dD/D[:,None]**2
    arg=s-q/D;clipped=(arg<p.lower)|(arg>p.sc)
    identity=np.zeros_like(dqs);identity[:,p.nx:]=np.eye(p.ns)
    dqs[clipped]=identity[clipped]
    return np.vstack([-H[:p.nx]*p.energy_scale/(p.length_scale*p.config['mechanical_equilibrium']['force_floor_N']),dqs])


def displacement_predictor(p,previous_x,previous_strain):
    x=np.array(previous_x,copy=True)
    x[:,2]+=(p.strain-previous_strain)*(p.family.X[:,2]-p.family.X[:,2].min())
    return p.prescribe(x)


def solve_stationary_branch_state(p,x_initial,s_initial):
    """Bounded SciPy residual solve; accepts stationary saddles as stationary."""
    c=settings();z0=p.pack(x_initial,s_initial);history=[];last=z0.copy()
    def fun(z):
        nonlocal last
        r=stationary_residual(p,z);last=z.copy()
        x,s=p.unpack(z);history.append(dict(evaluation=len(history),**p.metrics(x,s)))
        return r
    try:
        result=least_squares(fun,z0,jac=lambda z:residual_jacobian(p,z),
            bounds=(np.r_[np.full(p.nx,-np.inf),p.lower],np.r_[np.full(p.nx,np.inf),np.full(p.ns,p.sc)]),
            **c['least_squares'])
        z=result.x;success=bool(result.success);reason=str(result.message)
        solver=dict(status=int(result.status),nfev=int(result.nfev),njev=int(result.njev),cost=float(result.cost),optimality=float(result.optimality))
        exit_info=None
    except BranchDomainExit as exc:
        z=last;success=False;reason=str(exc);solver={};exit_info=dict(reason=reason,metrics=exc.metrics,z=exc.z.tolist())
    x,s=p.unpack(z);m=p.metrics(x,s)
    gates=dict(least_squares_converged=success,finite=bool(np.isfinite(z).all() and np.isfinite(m['total_J'])),
        force=m['force_residual']<=c['force_tolerance'],phase_KKT=m['projected_phase_KKT']<=c['phase_KKT_tolerance'],
        bounds=m['minimum_lower_margin']>=-c['bounds_tolerance'] and m['minimum_upper_margin']>0,
        irreversibility=bool(np.all(s>=p.lower)),high_damage=m['maximum_s']<p.sc,
        pure_tension=m['minimum_Green_strain']>=p.config['tolerances']['tensile_min_green_strain'],
        fixed_boundary=bool(np.array_equal(x[p.cloud.fixed],p.boundary[p.cloud.fixed])))
    return dict(x=x,s=s,z=z,z_initial=z0,lower=p.lower.copy(),accepted=all(gates.values()),gates=gates,
                reason=reason,solver=solver,exit=exit_info,metrics=m,history=history)


def active_sets(p,x,s):
    groups=p.classify(x,s,audit_config())
    result={'LOWER_STRONG_ACTIVE':groups['LOWER_BOUND_ACTIVE'], 'INTERIOR':groups['INTERIOR_ACTIVE'],
        'LOWER_WEAK_NEAR_SWITCH':groups['WEAK_LOWER_OR_UNCLASSIFIED'] & ~groups['NEAR_SC'],
        'NEAR_SC':groups['NEAR_SC']}
    return result,dict(counts={k:int(v.sum()) for k,v in result.items()},bond_IDs={k:np.flatnonzero(v).tolist() for k,v in result.items()},q_scaled=groups['q_scaled'].tolist())


def original_hessian_audit(label,p,x,s,eps):
    """Route only the old writer's destination; restore its global immediately.

    Used in independent processes, never threads. Numerical settings and the
    entire Hessian/audit implementation remain the previous verified code.
    """
    before=legacy_hessian.OUT
    try:
        legacy_hessian.OUT=OUT
        return legacy_hessian.hessian_audit((label,x,s,p.strain,p.lower,eps))
    finally:legacy_hessian.OUT=before


def probe_point(p,z,direction,amplitude):
    trial=z+amplitude*direction;x,s=p.unpack(trial)
    bound=bool(np.all(s>=p.lower) and np.all(s<p.sc))
    if not bound:return dict(amplitude=amplitude,feasible=False,reason='PHASE_BOUNDS')
    state=p.state(x,s);green=green_minimum(state)
    return dict(amplitude=amplitude,feasible=green>=p.config['tolerances']['tensile_min_green_strain'],
                minimum_Green_strain=green,reason=None if green>=p.config['tolerances']['tensile_min_green_strain'] else 'PURE_TENSION_DOMAIN')


def cone_audit(p,x,s,blocks):
    """Deterministic cone witnesses; never call positive sampling a cone proof."""
    c=settings();masks,_=active_sets(p,x,s);z=p.pack(x,s);_,grad=p.objective_gradient(z)
    H=blocks[1]['H_symmetric'];rng=np.random.default_rng(c['cone_seed']);candidates=[]
    strong=masks['LOWER_STRONG_ACTIVE'];weak=masks['LOWER_WEAK_NEAR_SWITCH'];interior=masks['INTERIOR']
    def add(name,v):
        norm=np.linalg.norm(v)
        if norm>1e-14:candidates.append((name,np.asarray(v)/norm))
    for sign in [-1,1]:add('reduced_minimum_eigenmode',sign*blocks[1]['minimum_mode_scaled'])
    for index in range(c['random_directions_per_category']):
        ux=rng.normal(size=p.nx);phase=rng.normal(size=p.ns)*interior
        lower_phase=np.abs(rng.normal(size=p.ns))*(strong|weak)
        add('displacement_only_random',np.r_[ux,np.zeros(p.ns)])
        add('phase_interior_random',np.r_[np.zeros(p.nx),phase])
        add('mixed_critical_random',np.r_[ux,phase+np.abs(rng.normal(size=p.ns))*weak])
        add('lower_bound_opening_auxiliary',np.r_[np.zeros(p.nx),lower_phase])
        add('mixed_lower_opening_auxiliary',np.r_[ux,phase+lower_phase])
    # Project negative eigenvectors onto weak-bound feasible directions, and use
    # the already-audited positive Hxx block to select their mechanical response.
    keep=np.r_[np.arange(p.nx),p.nx+np.flatnonzero(interior|weak)]
    vals,vecs=eigh(H[np.ix_(keep,keep)])
    for j in np.flatnonzero(vals<0)[:c['projected_negative_mode_limit']]:
        for sign in [-1,1]:
            v=np.zeros(len(z));v[keep]=sign*vecs[:,j]
            v[p.nx:][weak]=np.maximum(v[p.nx:][weak],0)
            if np.linalg.cond(H[:p.nx,:p.nx])<audit_config()['hxx_condition_limit']:
                v[:p.nx]=-solve(H[:p.nx,:p.nx],H[:p.nx,p.nx:]@v[p.nx:],assume_a='pos')
            add('projected_negative_mode_with_mechanical_response',v)
    rows=[];vectors=[]
    for index,(name,v) in enumerate(candidates):
        linear=float(grad@v);tangent=bool(np.all(abs(v[p.nx:][strong])<1e-14) and np.all(v[p.nx:][weak]>=-1e-14))
        critical=tangent and abs(linear)<=c['critical_linear_tolerance']
        probes=[probe_point(p,z,v,amp) for amp in c['feasibility_amplitudes']]
        rows.append(dict(index=index,category=name,linear_energy_variation=linear,critical_cone=critical,
            feasible=any(r['feasible'] for r in probes),curvatures=[float(v@b['H_symmetric']@v) for b in blocks],probes=probes))
        vectors.append(v)
    return rows,np.array(vectors),float(vals[0])


def audit_stationary_state(label,p,result):
    c=settings();x,s=result['x'],result['s'];rows=[];blocks=[]
    for eps in c['hessian_epsilon']:
        row=original_hessian_audit(label,p,x,s,eps);rows.append(row)
        b=read_npz(OUT/f'hessian/{label}/eps_{eps:g}/BLOCKS.npz');blocks.append(b)
    _,active=active_sets(p,x,s);cones,directions,full_critical_min=cone_audit(p,x,s,blocks)
    thresholds=np.array([r['curvature_threshold'] for r in rows]);eigs=np.array([r['minimum_reduced_eigenvalue'] for r in rows])
    signs=np.where(eigs>thresholds,1,np.where(eigs<-thresholds,-1,0))
    matrix_valid=all(r['Hessian_symmetry_relative_error']<=audit_config()['hessian_symmetry_relative_tolerance'] and
        r['Hss_analytic_relative_error']<=audit_config()['hss_relative_tolerance'] for r in rows)
    consistent=bool(np.all(signs==signs[0]) and matrix_valid)
    negative=[r for r in cones if r['critical_cone'] and r['feasible'] and np.all(np.array(r['curvatures'])<-thresholds)]
    if not result['accepted']:classification='UNRESOLVED'
    elif not consistent:classification='UNRESOLVED'
    elif negative:classification='UNSTABLE_SAMPLED_STATE'
    elif signs[0]==0:classification='MARGINAL'
    elif signs[0]>0 and (active['counts']['LOWER_WEAK_NEAR_SWITCH']==0 or full_critical_min>max(thresholds)):
        classification='STABLE_SAMPLED_STATE'
    else:classification='UNRESOLVED'
    selected=None;kind='reduced_eigenvector';v=blocks[1]['minimum_mode_scaled'];indicator=float(eigs[1])
    if negative:
        selected=min(negative,key=lambda r:np.mean(r['curvatures']))
        if signs[0]>0 or active['counts']['LOWER_WEAK_NEAR_SWITCH']:
            v=directions[selected['index']];kind='feasible_critical_cone_witness';indicator=float(selected['curvatures'][1])
    summary=dict(state=label,strain=p.strain,stationary=result['accepted'],classification=classification,
        status='STATIONARY_UNSTABLE' if classification=='UNSTABLE_SAMPLED_STATE' else classification,
        reduced_eigenvalues_min=eigs.tolist(),second_reduced_eigenvalues=[float(b['reduced_eigenvalues'][1]) for b in blocks],
        curvature_sign_resolved=consistent,curvature_sign_status='RESOLVED' if consistent else 'CURVATURE_SIGN_UNRESOLVED',
        stability_indicator=indicator,indicator_kind=kind,active=active,
        critical_cone_negative_witness=selected,full_critical_linear_space_minimum=full_critical_min,
        minimum_Hxx_eigenvalues=[r['minimum_mechanical_eigenvalue'] for r in rows],
        minimum_Schur_eigenvalues=[r['minimum_Schur_eigenvalue'] for r in rows],
        Hxx_condition_numbers=[r['Hxx_condition_number'] for r in rows],
        curvature_thresholds=thresholds.tolist(),metrics=result['metrics'],
        prescribed_displacement_m=p.strain*np.ptp(p.family.X[:,2]),
        reaction_N=float(-p.state(x,s).force[p.top,2].sum()))
    save_json(OUT/f'hessian/{label}/BRANCH_SUMMARY.json',summary)
    save_json(OUT/f'active_sets/{label}.json',active)
    save_json(OUT/f'hessian/{label}/CONE_PROBES.json',cones)
    np.savez_compressed(OUT/f'hessian/{label}/CONE_AND_MODE.npz',directions=directions,critical_mode=v,
                        x=x,s=s,lower=p.lower,scaled_state=p.pack(x,s))
    return summary


def align_modes(reference,current,nx):
    ref=np.asarray(reference);cur=np.asarray(current).copy()
    if ref@cur<0:cur=-cur
    def cos(a,b):
        den=np.linalg.norm(a)*np.linalg.norm(b)
        return None if den<1e-14 else float(a@b/den)
    return cur,dict(full=cos(ref,cur),displacement=cos(ref[:nx],cur[:nx]),phase=cos(ref[nx:],cur[nx:]))

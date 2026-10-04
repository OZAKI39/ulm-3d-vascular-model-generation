"""Diagnostic residual blocks and SciPy active-set solves; never a production solver."""
from scripts.c5_activation_audit import *
from scipy.optimize import least_squares, minimize
from scipy.sparse.linalg import LinearOperator, svds
from scipy.linalg import lu_factor, lu_solve, svdvals
import importlib


def original_once():
    folder=OUT/'original_solver';folder.mkdir(exist_ok=True)
    if (folder/'SUMMARY.json').exists():raise RuntimeError('Original solver is allowed exactly once')
    old=read(OLD/'branch_tracking/C5_P_005/STATIONARY.json')
    a=ca.read_npz(OLD/'branch_tracking/C5_P_005/STATE.npz');p=Problem(5,old['strain'],a['lower'])
    x,s=p.unpack(a['z_initial']);evaluations=[];termination=[];linear=[];capture={}
    dogbox=importlib.import_module('scipy.optimize._lsq.dogbox')
    lsq=tracker.least_squares;term=dogbox.check_termination;lsmr=dogbox.lsmr
    def logged_lsq(fun,z0,**kwargs):
        previous=[z0.copy()]
        def f(z):
            r=fun(z);xx,ss=p.unpack(z)
            evaluations.append(dict(evaluation=len(evaluations),x_norm=float(np.linalg.norm(z)),
                trial_step_norm=float(np.linalg.norm(z-previous[0])),mechanical_L2=float(np.linalg.norm(r[:p.nx])),
                mechanical_Linf=float(abs(r[:p.nx]).max()),phase_L2=float(np.linalg.norm(r[p.nx:])),
                phase_Linf=float(abs(r[p.nx:]).max()),cost=float(r@r/2),**p.metrics(xx,ss)))
            previous[0]=z.copy();return r
        result=lsq(f,z0,**kwargs);capture['result']=result;return result
    def logged_term(dF,F,dx,xnorm,ratio,ftol,xtol):
        status=term(dF,F,dx,xnorm,ratio,ftol,xtol)
        termination.append(dict(dF=float(dF),cost=float(F),step_norm=float(dx),x_norm=float(xnorm),ratio=float(ratio),
            ftol_threshold=float(ftol*F),xtol_threshold=float(xtol*(xtol+xnorm)),status=status,
            ftol_trigger=bool(dF<ftol*F and ratio>.25),xtol_trigger=bool(dx<xtol*(xtol+xnorm))))
        return status
    def logged_lsmr(A,b,**kwargs):
        r=lsmr(A,b,**kwargs)
        linear.append(dict(kwargs=kwargs,istop=int(r[1]),iterations=int(r[2]),normr=float(r[3]),normar=float(r[4]),
            norma=float(r[5]),conda_estimate=float(r[6]),normx=float(r[7])))
        return r
    start=time.time()
    with patch.object(tracker,'least_squares',logged_lsq),patch.object(dogbox,'check_termination',logged_term),patch.object(dogbox,'lsmr',logged_lsmr):
        answer=spectral.spectral_stationary(p,x,s)
    result=capture['result'];r=result.fun;J=result.jac;gradient=J.rmatvec(r)
    record={k:v for k,v in answer.items() if k not in ['x','s','z','z_initial','lower']}
    record.update(strain=p.strain,wall_seconds=time.time()-start,x_norm=float(np.linalg.norm(result.x)),
                  JT_r_L2=float(np.linalg.norm(gradient)),JT_r_Linf=float(abs(gradient).max()),
                  evaluations=evaluations,termination_checks=termination,lsmr_calls=linear,
                  instrumentation='Read-only wrappers around original least_squares, dogbox check_termination and lsmr; no option or residual changed.')
    save(folder/'SUMMARY.json',record)
    np.savez_compressed(folder/'STATE.npz',**{k:answer[k] for k in ['x','s','z','z_initial','lower']},residual=r,gradient=gradient)
    row,arrays=record_fields(p,answer['x'],answer['s']);save(folder/'KINK.json',row);np.savez_compressed(folder/'KINK.npz',**arrays)
    print('ORIGINAL',record['reason'],record['metrics']['force_residual'],record['metrics']['projected_phase_KKT'],flush=True)


class ResidualBlocks:
    """Compact Jacobian of force/q equations; not a stability Hessian.

    The 225 displacement columns are central differences of unchanged mechanics
    and Y. Phase-only and mixed derivatives follow the existing energy gradient.
    This stores O(Ns*Nx), never a dense Ns by Ns matrix.
    """
    def __init__(self,p,z,step=3e-6):
        self.p=p;self.z=np.array(z,copy=True);self.x,self.s=p.unpack(z)
        self.Y,self.q,self.D,self.qs=q_fields(p,self.x,self.s)
        nx=p.nx;ns=p.ns;floor=p.config['mechanical_equilibrium']['force_floor_N']
        self.M=np.empty((nx,nx));self.Yx=np.empty((ns,nx))
        for j in range(nx):
            dz=np.zeros_like(z);dz[j]=step
            xp,sp=p.unpack(z+dz);xm,sm=p.unpack(z-dz)
            plus=p.state(xp,sp);minus=p.state(xm,sm)
            self.M[:,j]=(plus.force[p.free]-minus.force[p.free]).ravel()/(2*step*floor)
            self.Yx[:,j]=(plus.W0-minus.W0)/(2*step)
        self.B=self.Yx.T*(2*(1-self.s)*p.weights/(p.length_scale*floor))[None,:]
        self.Qx=-2*(1-self.s)[:,None]*self.Yx
        self.Qs=2*self.Y
        arg=self.s-self.qs;self.clipped=(arg<p.lower)|(arg>p.sc)
        self.C=self.Qx/self.D[:,None]-self.q[:,None]*self.Yx*(self.Y>p.yc)[:,None]/self.D[:,None]**2
        self.C[self.clipped]=0
        self.diag=np.where(self.clipped,1.,self.Qs/self.D)

    def operator(self):
        nx=self.p.nx;ns=self.p.ns
        def mv(v):
            v=np.asarray(v).ravel();u,s=v[:nx],v[nx:]
            return np.r_[self.M@u+self.B@s,self.C@u+self.diag*s]
        def rm(v):
            v=np.asarray(v).ravel();u,s=v[:nx],v[nx:]
            return np.r_[self.M.T@u+self.C.T@s,self.B.T@u+self.diag*s]
        return LinearOperator((nx+ns,nx+ns),matvec=mv,rmatvec=rm,dtype=float)

    def inverse_operator(self):
        # Linear diagnostic factorization of the 225x225 Schur block only.
        # It is never used to implement a nonlinear Newton step.
        p=self.p;C_over_d=self.C/self.diag[:,None]
        schur=self.M-self.B@C_over_d;lu=lu_factor(schur)
        def mv(v):
            v=np.asarray(v).ravel();a,b=v[:p.nx],v[p.nx:]
            x=lu_solve(lu,a-self.B@(b/self.diag));s=(b-self.C@x)/self.diag
            return np.r_[x,s]
        def rm(v):
            v=np.asarray(v).ravel();a,b=v[:p.nx],v[p.nx:]
            x=lu_solve(lu,a-self.C.T@(b/self.diag),trans=1);s=(b-self.B.T@x)/self.diag
            return np.r_[x,s]
        return LinearOperator((p.nx+p.ns,)*2,matvec=mv,rmatvec=rm,dtype=float)

    def reduced(self,active):
        a=np.asarray(active,int)
        return np.block([[self.M,self.B[:,a]],
                         [self.Qx[a]/self.p.yc,np.diag(self.Qs[a]/self.p.yc)]])


def validate_blocks(p,z,blocks,label):
    cfg=read(OLD/'provenance/PREDECLARED_CONFIG.json')
    with patch.object(tracker,'settings',lambda:cfg):oldJ=v2.matrix_free_residual_jacobian(p,z)
    J=blocks.operator();rng=np.random.default_rng(5910);rows=[]
    for i in range(4):
        v=rng.normal(size=len(z));v/=np.linalg.norm(v);w=rng.normal(size=len(z));w/=np.linalg.norm(w)
        # At lower-bound components, use a positive direction and one-sided
        # second-order probes. Preserve the selected clip branch in this check.
        v[p.nx:]=abs(v[p.nx:]);v/=np.linalg.norm(v)
        a=J@v;b=oldJ@v
        rows.append(dict(direction=i,old_Jv_relative_error=float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30)),
            adjoint_error=float(abs(w@(J@v)-v@J.rmatvec(w))/max(np.linalg.norm(a),1e-30))))
    result=dict(label=label,rows=rows,passed=max(r['old_Jv_relative_error'] for r in rows)<5e-5)
    save(OUT/f'jacobian/{label}/VALIDATION.json',result)
    assert result['passed'],'RESIDUAL_JACOBIAN_ADAPTER_VALIDATION_FAILED'
    return result


def conditioning(p,z,label):
    folder=OUT/f'jacobian/{label}';start=time.time();blocks=ResidualBlocks(p,z);J=blocks.operator()
    validate_blocks(p,z,blocks,label)
    r=tracker.stationary_residual(p,z,check_domain=False);rx=np.r_[r[:p.nx],np.zeros(p.ns)];rs=np.r_[np.zeros(p.nx),r[p.nx:]]
    rownorm=np.sqrt(np.r_[np.sum(blocks.M**2,axis=1)+np.sum(blocks.B**2,axis=1),np.sum(blocks.C**2,axis=1)+blocks.diag**2])
    colnorm=np.sqrt(np.r_[np.sum(blocks.M**2,axis=0)+np.sum(blocks.C**2,axis=0),np.sum(blocks.B**2,axis=0)+blocks.diag**2])
    mechanical=svdvals(blocks.M)
    result=dict(label=label,strain=p.strain,mechanical_sigma_max=float(mechanical[0]),mechanical_sigma_min=float(mechanical[-1]),
        mechanical_condition=float(mechanical[0]/mechanical[-1]),phase_diagonal_min=float(blocks.diag.min()),phase_diagonal_max=float(blocks.diag.max()),
        cross_xs_Frobenius=float(np.linalg.norm(blocks.B)),cross_sx_Frobenius=float(np.linalg.norm(blocks.C)),
        mechanical_residual_L2=float(np.linalg.norm(rx)),phase_residual_L2=float(np.linalg.norm(rs)),
        mechanical_gradient_L2=float(np.linalg.norm(J.rmatvec(rx))),phase_gradient_L2=float(np.linalg.norm(J.rmatvec(rs))),
        total_gradient_L2=float(np.linalg.norm(J.rmatvec(r))),
        row_norm_percentiles=np.percentile(rownorm,[0,5,50,95,100]).tolist(),column_norm_percentiles=np.percentile(colnorm,[0,5,50,95,100]).tolist(),
        storage_bytes=int(sum(a.nbytes for a in [blocks.M,blocks.B,blocks.C,blocks.Yx])),
        method='SciPy svds LM on compact projected residual J; smallest via LM on exact block inverse; validate both singular triplet residuals on original J.')
    try:
        u,high,vh=svds(J,k=1,which='LM',tol=1e-9,maxiter=500,rng=np.random.default_rng(21))
        inv=blocks.inverse_operator();ui,hi,vhi=svds(inv,k=1,which='LM',tol=1e-9,maxiter=500,rng=np.random.default_rng(22))
        lo=1/hi[0];ul=vhi[0];vl=ui[:,0]
        higherr=max(np.linalg.norm(J@vh[0]-high[0]*u[:,0]),np.linalg.norm(J.rmatvec(u[:,0])-high[0]*vh[0]))/high[0]
        lowerr=max(np.linalg.norm(J@vl-lo*ul),np.linalg.norm(J.rmatvec(ul)-lo*vl))/lo
        resolved=higherr<=1e-6 and lowerr<=1e-6
        result.update(sigma_max=float(high[0]),sigma_min=float(lo),sigma_max_relative_residual=float(higherr),sigma_min_relative_residual=float(lowerr),
                      condition_estimate=float(high[0]/lo) if resolved else None,condition_resolved=bool(resolved))
    except Exception as exc:
        result.update(condition_estimate=None,condition_resolved=False,error=str(exc))
    result['wall_seconds']=time.time()-start
    save(folder/'SUMMARY.json',result)
    np.savez_compressed(folder/'BLOCKS.npz',M=blocks.M,B=blocks.B,C=blocks.C,diag=blocks.diag,Yx=blocks.Yx,
                        row_norms=rownorm,column_norms=colnorm,z=z,gradient_mechanical=J.rmatvec(rx),gradient_phase=J.rmatvec(rs))
    print('CONDITION',label,result.get('condition_estimate'),result.get('sigma_min_relative_residual'),flush=True)
    return result


def jacobian_audit():
    a=ca.read_npz(OUT/'reproduction/C5_W_004.npz');r=read(OUT/'reproduction/C5_W_004.json');p=Problem(5,r['strain'],a['lower'])
    rows=[conditioning(p,p.pack(a['x'],a['s']),'A_last_accepted')]
    r=read(OUT/'activation/POINT.json');a=ca.read_npz(OUT/f'intact/{r["label"]}/STATE.npz');p=Problem(5,r['strain'])
    rows.append(conditioning(p,p.pack(a['x'],a['s']),'B_first_weak'))
    r=read(OUT/'original_solver/SUMMARY.json');a=ca.read_npz(OUT/'original_solver/STATE.npz');p=Problem(5,r['strain'],a['lower'])
    rows.append(conditioning(p,a['z'],'C_failed_coupled'))
    save(OUT/'jacobian/ALL.json',rows)


def selections(p,x,s):
    _,q,_,qs=q_fields(p,x,s)
    a1=np.flatnonzero(q<0);a2=np.flatnonzero(qs<=1e-7)
    union=np.union1d(a1,a2);nodes=np.unique(p.family.pairs[union])
    adjacent=np.any(np.isin(p.family.pairs,nodes),axis=1)&(qs<=1e-3)
    return dict(A1=a1,A2=a2,A3=np.union1d(union,np.flatnonzero(adjacent)))


class ReducedSystem:
    def __init__(self,p,active):
        self.p=p;self.active=np.asarray(sorted(set(active)),int)
    def expand(self,u):
        z=np.r_[u[:self.p.nx],self.p.lower.copy()];z[self.p.nx+self.active]=u[self.p.nx:];return z
    def compress(self,z):return np.r_[z[:self.p.nx],z[self.p.nx+self.active]]
    def residual(self,u):
        p=self.p;z=self.expand(u);x,s=p.unpack(z);state=p.state(x,s);_,q,_,_=q_fields(p,x,s)
        return np.r_[state.force[p.free].ravel()/p.config['mechanical_equilibrium']['force_floor_N'],q[self.active]/p.yc]
    def jacobian(self,u):return ResidualBlocks(self.p,self.expand(u)).reduced(self.active)


def full_acceptance(p,z,active):
    x,s=p.unpack(z);Y,q,D,qs=q_fields(p,x,s);metrics=p.metrics(x,s);inactive=np.ones(p.ns,bool);inactive[active]=False
    view=adr.boundary_view(p.cloud,p.top,p.strain);safe=True
    try:check_positions(view,x,GUARD)
    except Exception:safe=False
    state=p.state(x,s);J=np.linalg.det(state.Ftilde)
    active_res=float(np.max(abs(qs[active]))) if len(active) else 0.
    inactive_violation=float(np.max(np.maximum(-qs[inactive],0))) if inactive.any() else 0.
    gates=dict(force=metrics['force_residual']<=1e-9,KKT=metrics['projected_phase_KKT']<=1e-8,
        active_equalities=active_res<=1e-8,inactive_inequalities=inactive_violation<=1e-8,
        bounds=bool(np.all(s>=0) and np.all(s<p.sc)),irreversibility=bool(np.all(s>=p.lower)),
        inactive_exact=bool(np.array_equal(s[inactive],p.lower[inactive])),
        fixed_boundary=bool(np.array_equal(x[p.cloud.fixed],p.boundary[p.cloud.fixed])),
        finite_J=bool(np.isfinite(J).all() and np.all(J>0)),safety=safe)
    return dict(valid=bool(all(gates.values())),positive_phase=bool(s.max()>1e-8),gates=gates,metrics=metrics,
                active_scaled_q_max=active_res,inactive_scaled_violation_max=inactive_violation,
                active_size=len(active),positive_phase_count=int((s>1e-8).sum()),reaction_stress_Pa=float(-state.force[p.top,2].sum()/v2.L**2))


def update_active(p,z,active):
    x,s=p.unpack(z);_,_,_,qs=q_fields(p,x,s)
    inactive=np.ones(p.ns,bool);inactive[active]=False
    add=np.flatnonzero(inactive&(qs<-1e-8))
    remove=np.asarray([b for b in active if s[b]<=p.lower[b]+1e-8 and qs[b]>1e-8],int)
    return np.setdiff1d(np.union1d(active,add),remove),add,remove


def reduced_solve(p,x,s,active,label):
    folder=OUT/f'reduced/{label}';folder.mkdir(parents=True,exist_ok=True);start=time.time();history=[];z=p.pack(x,s)
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json')['reduced_solver'];changes=[]
    for iteration in range(11):
        system=ReducedSystem(p,active);u0=system.compress(z)
        u0[p.nx:]=np.maximum(u0[p.nx:],p.lower[system.active])
        def fun(u):
            r=system.residual(u)
            history.append(dict(active_set_iteration=iteration,evaluation=len(history),force=float(np.linalg.norm(r[:p.nx])),
                                active_q_Linf=float(abs(r[p.nx:]).max()) if len(active) else 0.,cost=float(r@r/2)))
            return r
        result=least_squares(fun,u0,jac=system.jacobian,
            bounds=(np.r_[np.full(p.nx,-np.inf),p.lower[system.active]],np.r_[np.full(p.nx,np.inf),np.full(len(active),p.sc)]),**cfg)
        z=system.expand(result.x);check=full_acceptance(p,z,system.active)
        new,add,remove=update_active(p,z,system.active)
        changes.append(dict(iteration=iteration,active_IDs=system.active.tolist(),added=add.tolist(),removed=remove.tolist(),
            solver_status=int(result.status),solver_message=str(result.message),nfev=int(result.nfev),njev=int(result.njev),
            cost=float(result.cost),optimality=float(result.optimality),**check))
        save(folder/'ITERATIONS.json',changes);save(folder/'HISTORY.json',history)
        print('REDUCED',label,iteration,'A',len(active),'valid',check['valid'],'force',check['metrics']['force_residual'],'KKT',check['metrics']['projected_phase_KKT'],flush=True)
        if check['valid'] or np.array_equal(new,system.active) or iteration==10:break
        active=new
    x,s=p.unpack(z);answer=dict(label=label,strain=p.strain,wall_seconds=time.time()-start,**check,
                              active_IDs=system.active.tolist(),active_set_updates=len(changes)-1,solver_message=str(result.message))
    save(folder/'SUMMARY.json',answer);np.savez_compressed(folder/'STATE.npz',x=x,s=s,z=z,lower=p.lower,active=system.active)
    return answer,z


def reduced_jacobian_check(p,z,active,label):
    system=ReducedSystem(p,active);u=system.compress(z);J=system.jacobian(u);rng=np.random.default_rng(631);rows=[]
    for k in range(3):
        v=rng.normal(size=len(u));v[p.nx:]=abs(v[p.nx:]);v/=np.linalg.norm(v)
        r0=system.residual(u)
        for e in [1e-5,3e-6,1e-6]:
            fd=(-3*r0+4*system.residual(u+e*v)-system.residual(u+2*e*v))/(2*e)
            rows.append(dict(direction=k,epsilon=e,relative_error=float(np.linalg.norm(fd-J@v)/max(np.linalg.norm(J@v),1e-30))))
    record=dict(label=label,rows=rows,passed=max(r['relative_error'] for r in rows)<5e-5)
    save(OUT/f'reduced/{label}/JACOBIAN_CHECK.json',record)
    assert record['passed']


def roots():
    point=read(OUT/'activation/POINT.json')['strain'];row,a=intact(point+1e-6,'root_offset_1e-6');p=Problem(5,row['strain'])
    sets=selections(p,a['x'],a['s']);z=p.pack(a['x'],a['s'])
    save(OUT/'reduced/SELECTIONS.json',{k:v.tolist() for k,v in sets.items()})
    reduced_jacobian_check(p,z,sets['A1'],'first_A1')
    answers=[];valid=None
    for name,active in sets.items():
        # Use the existing AT1 law once as a bound-feasible initial guess only.
        seed=spectral.at1_constrained_update(p.lower,a['Y'],p.yc,p.sc).s_new
        seed[np.setdiff1d(np.arange(p.ns),active)]=0
        r,zr=reduced_solve(p,a['x'],seed,active,f'first_{name}');answers.append(r)
        if r['valid'] and r['positive_phase'] and valid is None:valid=(r,zr)
    save(OUT/'reduced/INITIAL_RESULTS.json',answers)
    if valid is None:
        save(OUT/'reduced/ROOT_SCAN.json',answers)
        print('REDUCED FIRST POINT NOT RESOLVED: fallback investigation required',flush=True);return
    first,zfirst=valid;prefix=[first];previous=zfirst;lower=p.unpack(zfirst)[1];oldstrain=p.strain
    for k in range(1,6):
        e=point+(k+1)*1e-6
        if (OUT/f'intact/root_offset_{k+1}e-6/SUMMARY.json').exists():rr,aa=intact(e,f'root_offset_{k+1}e-6')
        else:rr,aa=intact(e,f'continuation_intact_{k:02d}')
        pp=Problem(5,e,lower);xp=tracker.displacement_predictor(pp,p.unpack(previous)[0],oldstrain)
        active=np.union1d(np.flatnonzero(lower>0),selections(pp,aa['x'],lower)['A2'])
        r,znew=reduced_solve(pp,xp,lower,active,f'continuation_{k:02d}');prefix.append(r)
        if not (r['valid'] and r['positive_phase']):break
        p=pp;previous=znew;lower=p.unpack(znew)[1];oldstrain=e
    save(OUT/'reduced/CONTINUATION.json',prefix)
    # Additional independent nearby existence queries keep the original s=0
    # previous-load bound. They are not additional continuation steps.
    scan=[first]
    for offset in [2e-6,5e-6,1e-5,2e-5,5e-5]:
        e=point+offset
        existing=next((r for r in prefix if abs(r['strain']-e)<1e-14 and r['valid']),None)
        if existing:
            scan.append(dict(existing,reused_continuation=True));continue
        rr,aa=intact(e,f'root_scan_{offset:g}');pp=Problem(5,e)
        active=selections(pp,aa['x'],aa['s'])['A2'];seed=spectral.at1_constrained_update(pp.lower,aa['Y'],pp.yc,pp.sc).s_new
        r,_=reduced_solve(pp,aa['x'],seed,active,f'offset_{offset:g}');scan.append(r)
    save(OUT/'reduced/ROOT_SCAN.json',scan)
    # Same-load comparison is required: no claim based on changing load step.
    e=read(OUT/'original_solver/SUMMARY.json')['strain'];rr,aa=intact(e,'original_failure_load');pp=Problem(5,e)
    active=selections(pp,aa['x'],aa['s'])['A2'];seed=spectral.at1_constrained_update(pp.lower,aa['Y'],pp.yc,pp.sc).s_new
    r,_=reduced_solve(pp,aa['x'],seed,active,'same_load_as_original')
    save(OUT/'reduced/SAME_LOAD.json',r)


def diagnostic_oracles():
    """Replay one frozen linear step and one optional monolithic energy solve."""
    source=read(OUT/'original_solver/SUMMARY.json');a=ca.read_npz(OUT/'original_solver/STATE.npz')
    p=Problem(5,source['strain'],a['lower']);z=a['z'];dogbox=importlib.import_module('scipy.optimize._lsq.dogbox')
    cfg=read(OLD/'provenance/PREDECLARED_CONFIG.json')
    with patch.object(tracker,'settings',lambda:cfg):J=v2.matrix_free_residual_jacobian(p,z)
    r=tracker.stationary_residual(p,z,check_domain=False);g=J.rmatvec(r)
    lb=np.r_[np.full(p.nx,-np.inf),p.lower];ub=np.r_[np.full(p.nx,np.inf),np.full(p.ns,p.sc)]
    on=np.zeros(len(z));on[z==lb]=-1;on[z==ub]=1
    blocked=on*g<0;free=~blocked;gcopy=g.copy();gcopy[blocked]=0
    op=dogbox.lsmr_operator(J,np.ones(len(z)),blocked)
    ls=dogbox.lsmr(op,r);newton=-ls[0][free]
    qa,qb=dogbox.build_quadratic_1d(J,gcopy,-gcopy)
    # Recorded steps never hit the original trust boundary, so Delta retains
    # its initialization norm. This is one frozen linear diagnostic, no new
    # nonlinear trial, continuation, or rerun of the original solver.
    delta=float(abs(a['z_initial']).max())
    step,bounds,hit=dogbox.dogleg_step(z[free],newton,gcopy[free],qa,qb,np.full(free.sum(),delta),lb[free],ub[free])
    negative=np.flatnonzero(free)[(z[free]==lb[free])&(newton<0)]
    fullstep=np.zeros(len(z));fullstep[free]=step
    replay=dict(replayed_step_norm=float(np.linalg.norm(step)),original_final_step_norm=source['termination_checks'][-1]['step_norm'],
        newton_step_norm=float(np.linalg.norm(newton)),negative_newton_step_at_lower_bound_count=len(negative),
        negative_newton_step_at_lower_bound_IDs=negative.tolist(),
        negative_phase_bond_IDs=(negative[negative>=p.nx]-p.nx).tolist(),
        original_active_bound_count=int(blocked.sum()),trust_bound_hit=bool(hit),
        lsmr_linear_relative_residual=float(np.linalg.norm(op@ls[0]-r)/np.linalg.norm(r)),
        lsmr_istop=int(ls[1]),lsmr_iterations=int(ls[2]),scope='Frozen final iterate, one linear solve plus existing SciPy dogleg_step; no original nonlinear solver rerun.')
    save(OUT/'original_solver/FROZEN_LINEAR_STEP_REPLAY.json',replay)
    print('LINEAR STEP REPLAY',replay['replayed_step_norm'],'negative lower components',len(negative),flush=True)
    same=read(OUT/'reduced/SAME_LOAD.json')
    if not same['valid']:
        save(OUT/'oracles/ENERGY.json',dict(status='NOT_NEEDED',reason='No valid same-load reduced root contradicting the original solve'))
        return
    history=[];start=time.time()
    def fg(z):
        f,g=p.objective_gradient(z);x,s=p.unpack(z)
        history.append(dict(evaluation=len(history),objective=f,**p.metrics(x,s)))
        return f,g
    result=minimize(fg,a['z_initial'],jac=True,method='L-BFGS-B',bounds=list(zip(lb,ub)),
        options=dict(ftol=1e-15,gtol=1e-13,maxiter=400,maxls=30,maxcor=20))
    x,s=p.unpack(result.x);active=np.flatnonzero(s>0);check=full_acceptance(p,result.x,active)
    ref=ca.read_npz(OUT/'reduced/same_load_as_original/STATE.npz')
    save(OUT/'oracles/ENERGY.json',dict(status='VALID_STATIONARY_ROOT' if check['valid'] else 'NOT_PHYSICALLY_CONVERGED',
        optimizer_success=bool(result.success),message=str(result.message),nit=int(result.nit),nfev=int(result.nfev),
        wall_seconds=time.time()-start,**check,phase_max_difference=float(abs(s-ref['s']).max()),
        displacement_max_difference_m=float(abs(x-ref['x']).max()),energy_difference_J=float(p.metrics(x,s)['total_J']-same['metrics']['total_J']),
        curvature='NOT_ASSESSED; no claim of a certified local minimum',history=history))
    np.savez_compressed(OUT/'oracles/ENERGY_STATE.npz',x=x,s=s,z=result.x,lower=p.lower)
    print('ENERGY ORACLE',check['valid'],check['metrics']['force_residual'],check['metrics']['projected_phase_KKT'],flush=True)


def scaling_diagnostic():
    """Row equilibration of a frozen Jacobian only; physical residual unchanged."""
    rows=[]
    for source in read(OUT/'jacobian/ALL.json'):
        name=source['label'];a=ca.read_npz(OUT/f'jacobian/{name}/BLOCKS.npz');p=Problem(5,source['strain'])
        b=ResidualBlocks.__new__(ResidualBlocks);b.p=p
        for key in ['M','B','C','diag']:setattr(b,key,a[key])
        J=b.operator();inv=b.inverse_operator();norm=a['row_norms'];shape=J.shape
        R=LinearOperator(shape,matvec=lambda v:(J@np.asarray(v).ravel())/norm,
             rmatvec=lambda v:J.rmatvec(np.asarray(v).ravel()/norm),dtype=float)
        Ri=LinearOperator(shape,matvec=lambda v:inv@(norm*np.asarray(v).ravel()),
             rmatvec=lambda v:norm*inv.rmatvec(np.asarray(v).ravel()),dtype=float)
        u,hi,v=svds(R,k=1,which='LM',tol=1e-9,maxiter=500,rng=np.random.default_rng(67))
        ui,hii,vi=svds(Ri,k=1,which='LM',tol=1e-9,maxiter=500,rng=np.random.default_rng(68))
        low=1/hii[0];left=vi[0];right=ui[:,0]
        err=max(np.linalg.norm(R@right-low*left),np.linalg.norm(R.rmatvec(left)-low*right))/low
        rows.append(dict(label=name,original_condition=source['condition_estimate'],
            row_equilibrated_condition=float(hi[0]/low) if err<=1e-6 else None,smallest_triplet_relative_residual=float(err),
            scope='Frozen row-normalized Jacobian only. No rescaled root solve or formal residual change.'))
    save(OUT/'jacobian/SCALING_DIAGNOSTIC.json',rows)


if __name__=='__main__':
    original_once();jacobian_audit();roots();diagnostic_oracles();scaling_diagnostic()

"""Scoped reuse of the V2 tracker/stability helpers for the new constitutive model."""
from contextlib import ExitStack,contextmanager
from scripts.spectral_split_audit import *
from scripts import spatial_stability_v2 as audit
from scipy.sparse.linalg import LinearOperator,eigsh,cg,minres
from scipy.linalg import solve as dense_solve

class HessianBlockView:
    """Expose existing V2 HVP through the block interface expected by its cone audit."""
    __array_priority__=10000
    def __init__(self,p,z,epsilon):
        self.p=p;self.op=v2.hessian_operator(p,z,epsilon);nx=p.nx;ns=p.ns
        self.xx=LinearOperator((nx,nx),matvec=lambda u:(self.op@np.r_[np.asarray(u).ravel(),np.zeros(ns)])[:nx],dtype=float)
        self.xs=LinearOperator((nx,ns),matvec=lambda q:(self.op@np.r_[np.zeros(nx),np.asarray(q).ravel()])[:nx],
            rmatvec=lambda u:(self.op@np.r_[np.asarray(u).ravel(),np.zeros(ns)])[nx:],dtype=float)
    def __getitem__(self,key):
        a,b=key
        if a.start is None and a.stop==self.p.nx and b.start is None and b.stop==self.p.nx:return self.xx
        if a.start is None and a.stop==self.p.nx and b.start==self.p.nx and b.stop is None:return self.xs
        raise TypeError('This audit requests only mechanical and cross blocks')
    def __rmatmul__(self,v):return self.op.rmatvec(v)

def condensed_solve(A,b,assume_a=None):
    """Dispatch existing cone condensation to mature SciPy CG for an operator."""
    if not isinstance(A,LinearOperator):return dense_solve(A,b,assume_a=assume_a)
    b=np.asarray(b);bn=np.linalg.norm(b)
    if bn==0:return np.zeros_like(b)
    x,info=cg(A,b/bn,rtol=1e-9,atol=0.,maxiter=10*A.shape[0]);x*=bn
    residual=float(np.linalg.norm(A@x-b)/bn)
    if info or residual>1e-7:raise RuntimeError(f'Mechanical condensation unresolved: info={info}, residual={residual}')
    return x

def fine_stability(p,result,label):
    """Reuse V2 cone candidates/feasibility/curvature with operator block storage."""
    assert result['accepted']
    pp=SpectralProblem(p.n,p.strain,p.lower,common_scale=True);z=pp.pack(result['x'],result['s']);masks,active=tracker.active_sets(pp,result['x'],result['s'])
    blocks=[];rows=[];rng=np.random.default_rng(195);mech_v0=rng.normal(size=pp.nx)
    for eps in v2.EPS:
        H=HessianBlockView(pp,z,eps)
        low,vec=eigsh(H.xx,k=1,which='SA',tol=1e-7,v0=mech_v0)
        high=eigsh(H.op,k=1,which='LM',tol=1e-7,v0=np.random.default_rng(193).normal(size=len(z)),return_eigenvectors=False)
        cond=float(abs(high[0])/low[0]);row=dict(epsilon=eps,Hxx_min=float(low[0]),coupled_spectral_radius_estimate=float(abs(high[0])),
            Hxx_condition_upper_estimate=cond,curvature_threshold=1e-7*max(abs(float(high[0])),1.),
            mechanical_eigen_residual=float(np.linalg.norm(H.xx@vec[:,0]-low[0]*vec[:,0])))
        blocks.append(H);rows.append(row)
    valid=all(r['Hxx_min']>r['curvature_threshold'] and r['Hxx_condition_upper_estimate']<1e10 and r['mechanical_eigen_residual']<1e-6 for r in rows)
    tests=[];vectors=[];selected=None;classification='STATIONARY_UNRESOLVED';localization=None
    if valid and active['counts']['LOWER_STRONG_ACTIVE']==pp.ns:
        classification='STATIONARY_STABLE'
    elif valid:
        original_directions=audit.phase_directions
        def directions(problem,H=None,op=None):return original_directions(problem,H=None,op=H.op if isinstance(H,HessianBlockView) else op)
        with patch.object(audit,'phase_directions',directions),patch.object(audit,'solve',condensed_solve),patch.object(audit,'affine',lambda problem:result['x']):
            tests,vectors,_=audit.cone_probes(pp,z,blocks,rows)
        negative=[r for r in tests if r['verified_negative']]
        if negative:
            selected=min(negative,key=lambda r:r['curvatures'][1]);classification='STATIONARY_UNSTABLE'
            v=vectors[selected['index']]
            # V2 localization is reused; provide actual nonaffine state to its activity diagnostic.
            original_active=tracker.active_sets
            with patch.object(tracker,'active_sets',lambda problem,x,s:original_active(problem,result['x'],result['s'])):
                localization=audit.mode_data(pp,v,label)
            np.savez_compressed(OUT/f'modes/C{p.n}_SELECTED.npz',phase=v[pp.nx:],midpoints=pp.family.X[pp.family.pairs].mean(axis=1),vector=v)
        else:
            keep=np.r_[np.arange(pp.nx),pp.nx+np.flatnonzero(~masks['LOWER_STRONG_ACTIVE'])]
            def reduced(H):
                def mv(w):
                    v=np.zeros(len(z));v[keep]=np.asarray(w).ravel();return (H.op@v)[keep]
                return LinearOperator((len(keep),len(keep)),matvec=mv,dtype=float)
            vals=[]
            for H,row in zip(blocks,rows):
                ev=eigsh(reduced(H),k=1,which='SA',tol=1e-7,v0=np.random.default_rng(149).normal(size=len(keep)),return_eigenvectors=False)
                vals.append(float(ev[0]));row['full_critical_space_minimum']=float(ev[0])
            if all(v>r['curvature_threshold'] for v,r in zip(vals,rows)):classification='STATIONARY_STABLE'
    save(OUT/f'cone_stability/{label}.json',dict(classification=classification,epsilons=rows,candidates=tests,
        storage='Existing V2 bound-aware HVP LinearOperator; no dense coupled Hessian',critical_cone_implementation='Unmodified spatial_stability_v2.cone_probes'))
    return dict(status=classification,negative_witness=selected,localization=localization)

def fine_state(n,strain,lower,x_initial,s_initial,label,do_stability=False):
    p=SpectralProblem(n,strain,lower);start=time.time();result=spectral_stationary(p,x_initial,s_initial)
    folder=OUT/f'branch_tracking/{label}';folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/'STATE.npz',**{k:result[k] for k in ['x','s','z','z_initial','lower']})
    record={k:v for k,v in result.items() if k not in ['x','s','z','z_initial','lower']};record.update(level=f'C{n}',strain=strain,wall_seconds=time.time()-start)
    save(folder/'STATIONARY.json',record)
    _,active=tracker.active_sets(p,result['x'],result['s']);active['fractions']={k:v/p.ns for k,v in active['counts'].items()}
    active['bond_midpoints_m']=p.family.X[p.family.pairs].mean(axis=1).tolist();save(OUT/f'active_sets/{label}.json',active)
    state=p.state(result['x'],result['s']);R=float(-state.force[p.top,2].sum())
    row=dict(label=label,level=f'C{n}',strain=strain,stationary=result['accepted'],gates=result['gates'],metrics=result['metrics'],
        reason=result['reason'],active_counts=active['counts'],active_fractions=active['fractions'],reaction_N=R,nominal_stress_Pa=R/v2.L**2,
        status='STATIONARY_UNRESOLVED',negative_witness=None,localization=None,stability_audited=False)
    if result['accepted'] and do_stability:
        row.update(fine_stability(p,result,label),stability_audited=True)
    save(folder/'SUMMARY.json',row);print(label,'stationary',result['accepted'],row['status'],result['reason'],flush=True)
    return result,row

def frozen_phase_initializer(p,x,s,label):
    """Existing ADR at the unchanged previous-load phase, used only as an initial guess."""
    from pd_clot.pfpd_candidate import adr_adapter_at1 as adapter
    from pd_clot.regularization.equilibrium import ForceState,check_positions
    view=adapter.boundary_view(p.cloud,p.top,p.strain)
    guard=dict(safety=dict(maximum_displacement_m=.01,maximum_radius_ratio=1.),pipe=dict(radius_m=1.))
    def factory(family,material,phase):
        fixed_phase=np.array(phase,copy=True)
        def force(xx):
            check_positions(view,xx,guard);state=p.state(xx,fixed_phase);i,j=p.family.pairs.T
            eps=np.maximum(np.linalg.norm(xx[j]-xx[i],axis=1)/p.family.length-1,0.)
            return ForceState(state.force,np.zeros_like(xx),eps,state.energy,float(np.linalg.det(state.Ftilde).min()),
                np.zeros_like(xx),np.zeros(len(xx),dtype=bool))
        return force
    t=time.time()
    with patch.object(adapter,'FrozenPhaseForce',factory):
        result=adapter.equilibrate_tensile(p.cloud,p.family,p.material,p.top,x,s,p.strain,p.config)
    save(OUT/f'branch_tracking/{label}/FROZEN_PHASE_INITIALIZER.json',dict(converged=result.converged,summary=result.summary,
        wall_seconds=time.time()-t,role='Initial guess only; previous-load phase held fixed. Original branch solver and all acceptance gates still required.',
        phase_updates=0,solver_settings_changed=False))
    np.savez_compressed(OUT/f'branch_tracking/{label}/INITIALIZER.npz',x=result.x,s=s)
    return result.x

def fine_schedule(n,warm_start=False,phase_predictor=False):
    assert read(OUT/'GATES.json')['GATE_3']=='PASSED'
    data=ca.read_npz(OUT/f'equilibrium/C{n}/{STRAINS[0]:.8f}/STATE.npz');p=SpectralProblem(n,STRAINS[0]);x=data['x'];s=p.lower.copy()
    estimate=STRAINS[0]/np.sqrt(data['Psi_plus'].max()/p.yc)
    extended=sorted(set([STRAINS[0]]+[float(estimate*f) for f in [.5,.75,.9,.97,1.,1.03,1.06,1.12]]+WINDOW))
    prefix=f'C{n}_P' if phase_predictor else (f'C{n}_W' if warm_start else f'C{n}_A')
    save(OUT/f'branch_tracking/{prefix}_SCHEDULE.json',dict(original_window=WINDOW,actual_window=extended,
        low_strain_onset_estimate=estimate,extension_reason='Accepted intact equilibrium already has max Psi_plus/Yc >1 below 0.009; extend downward to track from inactive phase. Estimate guides sampling only, never defines an accepted threshold.'))
    rows=[];saved=[];oldstrain=STRAINS[0];previous_stable=None;first_unstable=None;barrier=None
    for index,e in enumerate(extended):
        p=SpectralProblem(n,e,s);xp=tracker.displacement_predictor(p,x,oldstrain)
        label=f'{prefix}_{index:03d}'
        previous=OUT/f'branch_tracking/C{n}_W_{index:03d}'
        reusable=phase_predictor and (previous/'STATIONARY.json').exists() and read(previous/'STATIONARY.json')['accepted']
        if reusable:
            result=read(previous/'STATIONARY.json');result.update(ca.read_npz(previous/'STATE.npz'));row=read(previous/'SUMMARY.json')
            assert np.array_equal(result['lower'],s) and abs(row['strain']-e)<1e-15
        else:
            cached_initializer=previous/'INITIALIZER.npz'
            if phase_predictor and cached_initializer.exists() and np.array_equal(ca.read_npz(cached_initializer)['s'],s):xp=ca.read_npz(cached_initializer)['x']
            elif warm_start and index:xp=frozen_phase_initializer(p,xp,s,label)
            seed=s
            if phase_predictor:
                seed=at1_constrained_update(s,p.state(xp,s).W0,p.yc,p.sc).s_new
                save(OUT/f'branch_tracking/{label}/PHASE_INITIALIZER.json',dict(method='Existing at1_constrained_update at fixed x, one call only',
                    maximum_phase_guess=float(seed.max()),previous_load_lower_bound_unchanged=True,solver_parameters_changed=False))
            result,row=fine_state(n,e,s,xp,seed,label)
        rows.append(row);saved.append(result)
        save(OUT/f'branch_tracking/C{n}_PATH.json',rows)
        if not result['accepted']:barrier=row;break
        active=row['active_counts'];needs_cone=active['LOWER_WEAK_NEAR_SWITCH']+active['INTERIOR']>0
        if needs_cone:
            # No coupled stability analysis was run before a legal weak/interior state existed.
            if previous_stable is None and len(rows)>1:
                rp=saved[-2];rr=rows[-2];lastp=SpectralProblem(n,rr['strain'],rp['lower'])
                rr.update(fine_stability(lastp,rp,rr['label']),stability_audited=True)
                save(OUT/f'branch_tracking/{rr["label"]}/SUMMARY.json',rr)
                if rr['status']=='STATIONARY_STABLE':previous_stable=(rp,rr)
                else:barrier=rr;break
            row.update(fine_stability(p,result,row['label']),stability_audited=True)
            save(OUT/f'branch_tracking/{row["label"]}/SUMMARY.json',row)
            if row['status']=='STATIONARY_UNSTABLE':first_unstable=(result,row);break
            if row['status']=='STATIONARY_UNRESOLVED':barrier=row;break
            previous_stable=(result,row)
        x=result['x'];s=result['s'];oldstrain=e
    if previous_stable and first_unstable and barrier is None:
        for k in range(12):
            lo,hi=previous_stable,first_unstable
            if hi[1]['strain']-lo[1]['strain']<=1e-5:break
            mid=(lo[1]['strain']+hi[1]['strain'])/2;p=SpectralProblem(n,mid,lo[0]['s'])
            xp=tracker.displacement_predictor(p,lo[0]['x'],lo[1]['strain'])
            label=f'{prefix}_refine_{k:02d}'
            if warm_start:xp=frozen_phase_initializer(p,xp,lo[0]['s'],label)
            seed=lo[0]['s']
            if phase_predictor:seed=at1_constrained_update(seed,p.state(xp,seed).W0,p.yc,p.sc).s_new
            result,row=fine_state(n,mid,lo[0]['s'],xp,seed,label,do_stability=True);rows.append(row)
            if row['status']=='STATIONARY_STABLE':previous_stable=(result,row)
            elif row['status']=='STATIONARY_UNSTABLE':first_unstable=(result,row)
            else:barrier=row;break
    interval=[previous_stable[1]['strain'],first_unstable[1]['strain']] if previous_stable and first_unstable and barrier is None else None
    answer=dict(level=f'C{n}',interval=interval,width=interval[1]-interval[0] if interval else None,
        complete=bool(interval and interval[1]-interval[0]<=1e-5),last_stable=previous_stable[1] if previous_stable else None,
        first_unstable=first_unstable[1] if first_unstable else None,unresolved_barrier=barrier,states=rows,
        reached_phase_neighborhood=any(r['stationary'] and r['active_counts']['INTERIOR']+r['active_counts']['LOWER_WEAK_NEAR_SWITCH']>0 for r in rows))
    save(OUT/f'critical_brackets/C{n}_A.json',answer);return answer

@contextmanager
def candidate_context():
    for name in ['critical_brackets','active_sets','matrix_free','cone_stability','modes','localization']:(OUT/name).mkdir(exist_ok=True)
    with ExitStack() as scope:
        for name,value in [('OUT',OUT),('Problem',SpectralProblem),('stationary',spectral_stationary)]:scope.enter_context(patch.object(audit,name,value))
        scope.enter_context(patch.object(tracker,'probe_point',mixed_probe_point))
        yield

def C3_regression():
    assert read(OUT/'GATES.json')['GATE_3']=='PASSED'
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json');cfg['schedules']={'A':WINDOW};save(OUT/'provenance/PREDECLARED_CONFIG.json',cfg)
    with candidate_context():result=audit.schedule(3,'A')
    old=read(v2.OUT/'critical_brackets/C3_A.json');comparisons=[]
    for key in ['last_stable','first_unstable']:
        if result[key] is None:continue
        new=result[key];previous=old[key];comparisons.append(dict(endpoint=key,strain_difference=new['strain']-previous['strain'],
            energy_relative_difference=abs(new['metrics']['total_J']/previous['metrics']['total_J']-1),
            reaction_relative_difference=abs(new['reaction_N']/previous['reaction_N']-1),
            negative_witness_verified=new['negative_witness']['verified_negative'] if new['negative_witness'] else None))
    passed=bool(result['complete'] and np.max(abs(np.array(result['interval'])-old['interval']))<1e-5 and
        all(r['energy_relative_difference']<1e-8 and r['reaction_relative_difference']<1e-8 for r in comparisons))
    save(OUT/'comparisons/C3_STABILITY_REGRESSION.json',dict(passed=passed,interval=result['interval'],previous_interval=old['interval'],comparisons=comparisons))
    print('C3 STABILITY REGRESSION',passed,result['interval'],flush=True)
    return passed

if __name__=='__main__':C3_regression()

"""C5 physical-energy stability audit; isolated writers and immutable inputs."""
from pathlib import Path
import ast,hashlib,inspect,json,os,time,resource,subprocess,multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.linalg import eigh,cho_factor,cho_solve
from scipy.optimize import minimize
from scripts import c5_activation_audit as old
from scripts import petsc_vi_audit as petsc
from scripts.run_monolithic_stability_audit import full_hessian
ROOT=old.ROOT;OUT=ROOT/'verification/pfpd_at1_c5_stability_v1';RESULTS=ROOT/'results/pfpd_at1_c5_stability_v1';CACHE=RESULTS/'cache'
read=old.read;save=old.save
EPS=[1e-5,3e-6,1e-6]


def settings():return read(OUT/'provenance/PREDECLARED_CONFIG.json')
def problem(e,lower):return old.spectral.SpectralProblem(5,e,lower,common_scale=True)
def checksum(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def declare():
    import psutil,platform,scipy
    c=read(old.v2.OUT/'provenance/PREDECLARED_CONFIG.json');h=old.ca.audit_config()
    mem=psutil.virtual_memory();physical=psutil.cpu_count(logical=False);reserve=.2*mem.total;parent=768*1024**2;worker=384*1024**2
    safe=max(1,int((mem.available-reserve-parent)//worker));workers=min(physical-2,safe)
    config=dict(task='PFPD_AT1_C5_STABILITY_V1',s_tol=c['active_s_tolerance'],q_tol=c['active_q_scaled_tolerance'],
        epsilons=EPS,raw_symmetry_tolerance=h['hessian_symmetry_relative_tolerance'],hss_tolerance=h['hss_relative_tolerance'],
        cross_symmetry_tolerance=h['hessian_symmetry_relative_tolerance'],FD_magnitude_relative_tolerance=c['hessian_epsilon_relative_tolerance'],
        curvature_relative_tolerance=c['curvature_relative_tolerance'],spectral_tol=1e-10,random_seed=24680,random_directions=10,
        direct_curvature_relative_tolerance=.02,direct_mode_amplitudes=[1e-6,3e-6,1e-5],
        spectral_smaller_amplitudes=[3e-7,1e-7],cone_max_starts=16,cone_seed=c['cone_seed'],max_critical_phase=512,
        reproduction_tolerances=dict(force_absolute=1e-12,KKT_absolute=1e-12,energy_relative=1e-9,reaction_relative=1e-8,phase_absolute=1e-12),
        scaling=dict(length_scale_m=c['stability_length_scale_m'],phase_scale=1.,energy_scale_J=c['stability_energy_scale_J'],source=str(old.v2.OUT/'provenance/PREDECLARED_CONFIG.json')),
        workers=workers,resource_budget=dict(physical_cores=physical,logical_cores=psutil.cpu_count(),RAM_total=mem.total,RAM_available=mem.available,reserve_RAM=reserve,parent_bytes=parent,worker_bytes=worker,memory_safe_workers=safe,BLAS_threads=1),
        continuation=dict(step=2.5e-5,max_strain=.006,max_s=.10,refinement_width=1e-6,max_refinements=10),
        curvature_definition='Dimensionless physical energy Hessian in one common scale. Threshold=1e-7*max(reduced Hessian spectral norm,1). Schur eigenvalue sign is a certificate, not full-mode Rayleigh magnitude.',
        direct_energy='Unmodified objective_gradient energy; raw energy differences and gradient-corrected estimates both saved. No energy formula copied.',
        negative_witness_rule='Feasible cone, negative Hessian curvatures at all three FD scales, direct physical-energy sign agreement at all valid nonswitching requested amplitudes; no negative search result is never stability proof.',
        implementation='Reuse exact existing full_hessian gradient-difference loop on reduced physical coordinates; process independent columns, no residual Jacobian used.',
        production_changes=False,GPU='NOT USED')
    save(OUT/'provenance/PREDECLARED_CONFIG.json',config);save(OUT/'SCALING.json',config['scaling'])
    save(OUT/'provenance/HARDWARE.json',dict(**config['resource_budget'],python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        cuda_backends={k:__import__('importlib.util',fromlist=['find_spec']).find_spec(k) is not None for k in ['cupy','jax','torch']},
        prior_PETSc=read(petsc.OUT/'provenance/PETSC_ENVIRONMENT.json')))
    print('DECLARED',workers,'workers',flush=True)


def saved_inputs():
    b=read(old.OUT/'activation/BRACKET.json');definitions=[]
    for label,name in [('S0_LAST',b['lower']['label']),('S_WEAK','activation_point')]:
        folder=old.OUT/f'intact/{name}';r=read(folder/'SUMMARY.json');definitions.append((label,r['strain'],folder/'STATE.npz',folder/'SUMMARY.json'))
    for i,t in enumerate(petsc.targets(),1):
        folder=petsc.OUT/f'continuation/PETSC_RSLS/{t["label"]}';definitions.append((f'S{i}',t['strain'],folder/'STATE.npz',folder/'SUMMARY.json'))
    return definitions


def classify(s,lower,qs,sc,cfg=None):
    cfg=cfg or settings();lo=abs(s-lower)<=cfg['s_tol'];up=abs(s-sc)<=cfg['s_tol'];qt=cfg['q_tol']
    return dict(STRONG_LOWER=lo&(qs>qt),WEAK_LOWER=lo&(abs(qs)<=qt),INTERIOR=(s>lower+cfg['s_tol'])&(s<sc-cfg['s_tol']),
        WEAK_UPPER=up&(abs(qs)<=qt),STRONG_UPPER=up&(qs<-qt))


def spectral_data(p,x,s):
    eigen=p.state(x,s).split_response.split.eigenvalues
    return dict(min_abs_principal_Green=float(abs(eigen).min()),near_zero_BAQP=int(np.any(abs(eigen)<settings()['spectral_tol'],axis=1).sum()),near_zero_eigenvalues=int((abs(eigen)<settings()['spectral_tol']).sum())),eigen


def reproduce():
    pre=read(OUT/'PRE_RUN_MANIFEST.json');rows=[];previous=None
    for label,e,path,summary_path in saved_inputs():
        for f in [path,summary_path]:
            assert checksum(f)==pre['files'][str(f.relative_to(ROOT))],'SAVED_STATIONARY_STATE_MISMATCH'
        a=old.ca.read_npz(path);ref=read(summary_path);lower=a.get('lower',np.zeros_like(a['s']));p=problem(e,lower);x,s=a['x'],a['s'];z=p.pack(x,s)
        result=petsc.acceptance(p,z,ref.get('reason',2));met=result['metrics'];oldmet=ref['metrics'];reaction=result['reaction_stress_Pa']
        errors=dict(force_absolute=abs(met['force_residual']-oldmet['force_residual']),KKT_absolute=abs(met['projected_phase_KKT']-oldmet['projected_phase_KKT']),
                    energy_relative=abs(met['total_J']/oldmet['total_J']-1),phase_absolute=abs(float(s.max())-oldmet['maximum_s']))
        if 'reaction_stress_Pa' in ref:errors['reaction_relative']=abs(reaction/ref['reaction_stress_Pa']-1)
        agreement=all(v<=settings()['reproduction_tolerances'][k] for k,v in errors.items())
        if label not in ['S0_LAST','S_WEAK','S1']:agreement=agreement and np.array_equal(lower,previous)
        if label=='S1':agreement=agreement and np.count_nonzero(lower)==0
        previous=s.copy();Y,q,D,qs=old.q_fields(p,x,s);groups=classify(s,lower,qs,p.sc);spec,_=spectral_data(p,x,s)
        passed=result['scientifically_accepted'] and agreement and np.all(np.sum(list(groups.values()),axis=0)==1)
        row=dict(label=label,strain=e,source=str(path),source_SHA256=checksum(path),summary_source=str(summary_path),summary_SHA256=checksum(summary_path),
                 passed=bool(passed),reproduction_errors=errors,reference_reaction='PREVIOUS_SUMMARY' if 'reaction_stress_Pa' in ref else 'NOT_PREVIOUSLY_STORED_IN_INTACT_SUMMARY',
                 **result,counts={k:int(v.sum()) for k,v in groups.items()},spectral=spec)
        save(OUT/f'states/{label}/INPUT_VERIFICATION.json',row);np.savez_compressed(OUT/f'states/{label}/STATE.npz',x=x,s=s,lower=lower,strain=e)
        rows.append(row);print('INPUT',label,passed,row['counts'],spec,flush=True)
        if not passed:save(OUT/'INPUT_VERIFICATION.json',dict(passed=False,states=rows));raise RuntimeError('SAVED_STATIONARY_STATE_MISMATCH')
    save(OUT/'INPUT_VERIFICATION.json',dict(passed=True,states=rows))


class CriticalProblem:
    """Restrict existing physical energy/gradient without reproducing its law."""
    def __init__(self,p,z,critical):
        self.p=p;self.base=z.copy();self.critical=np.asarray(critical,int);self.keep=np.r_[np.arange(p.nx),p.nx+self.critical]
        self.nx=p.nx;self.lower=p.lower[self.critical];self.sc=p.sc;self.probe_records=[]
        self.e0=p.state(*p.unpack(z)).split_response.split.eigenvalues.copy()
    def expand(self,z):
        full=self.base.copy();full[self.keep]=z;return full
    def objective_gradient(self,z):
        full=self.expand(z);value,g=self.p.objective_gradient(full)
        eigen=self.p.state(*self.p.unpack(full)).split_response.split.eigenvalues
        raw=np.sign(eigen)!=np.sign(self.e0)
        strong_switch=(self.e0*eigen<0)&(np.maximum(abs(self.e0),abs(eigen))>settings()['spectral_tol'])
        self.probe_records.append(dict(raw_sign_changes=int(raw.sum()),resolved_sign_changes=int(strong_switch.sum())))
        return value,g[self.keep]


def spectral_restriction(p,x,s,critical):
    e=p.state(x,s).split_response.split.eigenvalues
    near=np.any(abs(e)<settings()['spectral_tol'],axis=1)
    movable=np.zeros(p.ns,bool);movable[critical]=True
    exposed=near&((s!=0)|movable)
    return dict(near_zero_bonds=int(near.sum()),near_zero_with_unequal_split_weight=int((near&(s!=0)).sum()),
        near_zero_in_critical_phase=int((near&movable).sum()),exposed_switches=int(exposed.sum()),
        cancellation_identity='On intact strong-bound bonds ds=0 and g=1, so Psi_plus+Psi_minus is the original smooth SVK energy on the critical subspace. No split smoothing.',
        smooth_restricted_energy=bool(not exposed.any()),lame_lambda=float(p.material.lame))


def column_chunk(task):
    strain,lower,x,s,critical,eps,columns=task;p=problem(strain,lower);z=p.pack(x,s);proxy=CriticalProblem(p,z,critical);u=z[proxy.keep]
    # Reuse the inspected existing full_hessian loop verbatim, changing only
    # iteration indices to this process's independent columns.
    fn=ast.parse(inspect.getsource(full_hessian)).body[0];loop=next(n for n in fn.body if isinstance(n,ast.For));loop.iter=ast.Name(id='columns',ctx=ast.Load())
    n=len(u);H=np.zeros((n,n));g0=proxy.objective_gradient(u)[1];schemes=[]
    scope=dict(np=np,problem=proxy,z=u,epsilon=eps,H=H,g0=g0,schemes=schemes,columns=columns)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[loop],type_ignores=[])),'<unchanged physical full_hessian column loop>','exec'),scope)
    return columns,H[:,columns],schemes,proxy.probe_records,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024


def state_key(p,x,s):
    h=hashlib.sha256()
    for a in [x,s,p.lower]:h.update(np.asarray(a,dtype='<f8').tobytes())
    h.update(json.dumps(dict(strain=p.strain,physics=p.config,scaling=settings()['scaling'],kernel='physical_full_hessian_v1',eps=EPS),sort_keys=True).encode())
    h.update(inspect.getsource(full_hessian).encode());return h.hexdigest()


def build_blocks(p,x,s,critical,pool):
    key=state_key(p,x,s);folder=CACHE/key;folder.mkdir(exist_ok=True);z=p.pack(x,s);keep=np.r_[np.arange(p.nx),p.nx+critical];n=len(keep);rows=[];matrices=[]
    value,gradient=p.objective_gradient(z);np.savez_compressed(folder/'ENERGY_GRADIENT.npz',energy=value,gradient=gradient,z=z,critical=critical)
    analytic=2*p.weights[critical]*p.state(x,s).W0[critical]/p.energy_scale
    for eps in EPS:
        path=folder/f'BLOCKS_{eps:g}.npz';meta=folder/f'BLOCKS_{eps:g}.json'
        if path.exists() and meta.exists():
            record=read(meta);assert record['state_hash']==key;H=old.ca.read_npz(path)['H_raw'];record['cache_hit']=True
        else:
            start=time.perf_counter();chunks=[c.tolist() for c in np.array_split(np.arange(n),settings()['workers']) if len(c)]
            tasks=[(p.strain,p.lower,x,s,critical,eps,c) for c in chunks];H=np.empty((n,n));schemes=[None]*n;probes=[];peak=0
            for ids,a,schemes_chunk,probe,rss in pool.map(column_chunk,tasks):
                H[:,ids]=a
                for j,scheme in zip(ids,schemes_chunk):schemes[j]=scheme
                probes.extend(probe);peak=max(peak,rss)
            record=dict(state_hash=key,epsilon=eps,wall_seconds=time.perf_counter()-start,schemes=schemes,worker_peak_RSS_bytes=peak,
                raw_asymmetry=float(np.linalg.norm(H-H.T)/max(np.linalg.norm(H),1e-30)),Hxx_raw_asymmetry=float(np.linalg.norm(H[:p.nx,:p.nx]-H[:p.nx,:p.nx].T)/np.linalg.norm(H[:p.nx,:p.nx])),
                cross_relative_error=float(np.linalg.norm(H[:p.nx,p.nx:]-H[p.nx:,:p.nx].T)/max(np.linalg.norm(H[:p.nx,p.nx:]),1e-30)) if len(critical) else 0.,
                Hss_relative_error=float(np.linalg.norm(H[p.nx:,p.nx:]-np.diag(analytic))/max(np.linalg.norm(analytic),1e-30)) if len(critical) else 0.,
                raw_spectral_switches=sum(r['raw_sign_changes'] for r in probes),resolved_spectral_switches=sum(r['resolved_sign_changes'] for r in probes),cache_hit=False)
            np.savez_compressed(path,H_raw=H,Hss_analytic=analytic,critical=critical);save(meta,record)
        rows.append(record);matrices.append((H+H.T)/2)
    save(folder/'STATE_HASH.json',dict(hash=key,strain=p.strain,physics=p.config,scaling=settings()['scaling'],critical=critical.tolist()))
    return folder,matrices,rows


def spectral_switch(before,after):
    raw=np.sign(before)!=np.sign(after);tol=settings()['spectral_tol']
    return dict(raw_count=int(raw.sum()),resolved_count=int(((before*after<0)&(np.maximum(abs(before),abs(after))>tol)).sum()),
                status='SPECTRAL_SWITCH_DURING_PROBE' if raw.any() else 'NO_SPECTRAL_SWITCH')


def cone_feasible(pvec,groups,critical,tol=1e-12):
    return bool(np.all(pvec[groups['WEAK_LOWER'][critical]]>=-tol) and np.all(pvec[groups['WEAK_UPPER'][critical]]<=tol))


def direct_probe(p,z,v,curvature,amplitudes=None):
    from scripts.c5_precise_energy import energy,directional_value
    f64,g0=p.objective_gradient(z);x0,s0=p.unpack(z);f0=energy(p,x0,s0)
    energy_agreement=float(abs(f0/f64-1));assert energy_agreement<1e-12,"PHYSICAL_ENERGY_PRECISION_ADAPTER_MISMATCH"
    base=p.state(*p.unpack(z)).split_response.split.eigenvalues.copy();linear=float(g0@v);rows=[]
    for amp in amplitudes or settings()['direct_mode_amplitudes']:
        points=[]
        for sign in [-1,1]:
            t=sign*amp;trial=z+t*v;x,s=p.unpack(trial)
            feasible=bool(np.all(s>=p.lower) and np.all(s<p.sc))
            row=dict(sign=sign,amplitude=t,feasible=feasible)
            if feasible:
                value,gradient=p.objective_gradient(trial);e=p.state(x,s).split_response.split.eigenvalues
                raw_delta=value-f64;delta=directional_value(p,x0,s0,v,t)-f0;one=2*(delta-np.longdouble(t)*linear)/np.longdouble(amp)**2
                row.update(raw_float64_delta_J=float(raw_delta*p.energy_scale),delta_normalized_energy=float(delta),delta_J=float(delta*p.energy_scale),quadratic_reference_J=float(.5*amp**2*curvature*p.energy_scale),
                    linear_J=float(t*linear*p.energy_scale),one_sided_curvature=float(one),relative_error=float(abs(one-curvature)/max(abs(curvature),1e-30)),
                    spectral=spectral_switch(base,e),positive_J=bool(np.all(np.linalg.det(p.state(x,s).Ftilde)>0)))
            points.append(row)
        valid=[r for r in points if r['feasible']];both=len(valid)==2
        switch=any(r['spectral']['raw_count'] for r in valid)
        centered=(sum(r['delta_normalized_energy'] for r in valid)/amp**2) if both else None
        selected=[r['one_sided_curvature'] for r in valid] if switch or not both else [centered]
        rows.append(dict(amplitude=amp,points=points,centered_curvature=float(centered) if centered is not None else None,
            centered_used=both and not switch,validation_scheme='one_sided_physical_energy' if switch or not both else 'centered_physical_energy',
            selected_curvatures=selected,sign_consistent=bool(selected and all(c*curvature>0 for c in selected)),
            max_relative_error=max((abs(c-curvature)/max(abs(curvature),1e-30) for c in selected),default=None)))
    return dict(energy_oracle='Identical physical energy code with extended arithmetic',baseline_energy_relative_agreement=energy_agreement,linear_variation=linear,quadratic_curvature=float(curvature),rows=rows,
        passed=bool(all(r['sign_consistent'] for r in rows)),
        max_relative_error=max((r['max_relative_error'] for r in rows if r['max_relative_error'] is not None),default=None),
        rule='A switching amplitude is excluded from centered validation; admissible one-sided energy used. Raw sign changes include roundoff near exact zeros and are never silently discarded.')


def mechanical_oracle(p,z,H):
    rng=np.random.default_rng(settings()['random_seed']);rows=[];gradient_errors=[];f0,g=p.objective_gradient(z)
    for i in range(settings()['random_directions']):
        v=np.zeros_like(z);v[:p.nx]=rng.normal(size=p.nx);v/=np.linalg.norm(v);curvature=float(v[:p.nx]@H@v[:p.nx]);check=direct_probe(p,z,v,curvature,EPS)
        rows.append(dict(direction=i,**check))
        for eps in EPS:
            fd=(p.objective_gradient(z+eps*v)[0]-p.objective_gradient(z-eps*v)[0])/(2*eps)
            gradient_errors.append(float(abs(fd-g@v)/max(np.linalg.norm(g),1e-4)))
    return dict(passed=all(r['passed'] and r['max_relative_error']<=settings()['direct_curvature_relative_tolerance'] for r in rows),
        directions=rows,max_relative_error=max(r['max_relative_error'] for r in rows),gradient_force_identity=True,
        gradient_directional_max_error=max(gradient_errors))


def select_cone(S,groups,critical,growth):
    m=len(critical);vals,vec=eigh(S);minimum=vec[:,0];tol=1e-12;trials=[]
    def add(name,p):
        if np.linalg.norm(p)<=1e-16:return
        p=p/np.linalg.norm(p)
        if cone_feasible(p,groups,critical,tol):trials.append((float(p@S@p),name,p))
    add('unconstrained_minimum',minimum);add('unconstrained_minimum_flip',-minimum)
    def project(p):
        p=p.copy();lo=groups['WEAK_LOWER'][critical];up=groups['WEAK_UPPER'][critical]
        p[lo]=np.maximum(p[lo],0);p[up]=np.minimum(p[up],0);return p
    if vals[0]>0 and trials:return min(trials,key=lambda x:x[0]),dict(method='POSITIVE_UNCONSTRAINED_CERTIFICATE',lambda_min=float(vals[0]),optimizer_run=False)
    if trials and trials[0][0]<0:return min(trials,key=lambda x:x[0]),dict(method='FEASIBLE_UNCONSTRAINED_NEGATIVE',lambda_min=float(vals[0]),optimizer_run=False)
    starts=[project(minimum),project(-minimum),project(np.ones(m)),project(-np.ones(m)),project(growth)]
    if m<=8:starts.extend(project(q) for q in np.eye(m));starts.extend(project(-q) for q in np.eye(m))
    rng=np.random.default_rng(settings()['cone_seed'])
    while len(starts)<16:starts.append(project(rng.normal(size=m)))
    starts=[v/np.linalg.norm(v) for v in starts if np.linalg.norm(v)>1e-16][:16]
    # An SPD S already certifies the entire cone; use feasible probes only.
    for i,p in enumerate(starts):add(f'feasible_start_{i}',p)
    runs=[]
    if vals[0]<=0:
        bounds=[(0,None) if groups['WEAK_LOWER'][i] else (None,0) if groups['WEAK_UPPER'][i] else (None,None) for i in critical]
        for i,p in enumerate(starts):
            answer=minimize(lambda v:float(v@S@v),p,jac=lambda v:2*S@v,method='SLSQP',bounds=bounds,
                constraints=[dict(type='eq',fun=lambda v:float(v@v-1),jac=lambda v:2*v)],options=dict(ftol=1e-13,maxiter=200))
            v=answer.x/np.linalg.norm(answer.x);add(f'SLSQP_{i}',v);runs.append(dict(start=i,success=bool(answer.success),objective=float(v@S@v),constraint_error=float(abs(answer.x@answer.x-1))))
            if trials and min(r[0] for r in trials)<-1e-7:break
    return min(trials,key=lambda r:r[0]),dict(method='POSITIVE_UNCONSTRAINED_CERTIFICATE' if vals[0]>0 else 'CONSTRAINED_WITNESS_SEARCH',lambda_min=float(vals[0]),optimizer_run=bool(runs),runs=runs)


def audit_state(label,p,x,s,pool):
    start=time.perf_counter();cfg=settings();Y,q,D,qs=old.q_fields(p,x,s);groups=classify(s,p.lower,qs,p.sc)
    critical=np.flatnonzero(groups['WEAK_LOWER']|groups['INTERIOR']|groups['WEAK_UPPER']);m=len(critical);z=p.pack(x,s);spec,eigen=spectral_data(p,x,s);spec.update(spectral_restriction(p,x,s,critical))
    state_folder=OUT/f'states/{label}';state_folder.mkdir(parents=True,exist_ok=True)
    if m>cfg['max_critical_phase']:raise RuntimeError('LARGE_CRITICAL_PHASE_SET')
    folder,matrices,assembly=build_blocks(p,x,s,critical,pool);blocks=[];eigenvalues=[]
    # Raw reciprocity has been recorded before symmetric projection.
    for H,meta in zip(matrices,assembly):
        Hxx=H[:p.nx,:p.nx];Hx=H[:p.nx,p.nx:];Hss=np.diag(2*p.weights[critical]*Y[critical]/p.energy_scale)
        H[p.nx:,p.nx:]=Hss
        vx,ux=eigh(Hxx);norm=float(np.max(abs(eigh(H,eigvals_only=True))));threshold=cfg['curvature_relative_tolerance']*max(norm,1.)
        record=dict(epsilon=meta['epsilon'],lambda_x_min=float(vx[0]),lambda_x_max=float(vx[-1]),cond_Hxx=float(vx[-1]/vx[0]),curvature_tol=threshold,m=m)
        if vx[0]>threshold:
            t=time.perf_counter();factor=cho_factor(Hxx,lower=True);X=cho_solve(factor,Hx);S=Hss-Hx.T@X
            record.update(cholesky_and_solve_seconds=time.perf_counter()-t,S_bytes=S.nbytes,lambda_S_min=float(eigh(S,eigvals_only=True)[0]) if m else None,
                condensation_residual=float(np.linalg.norm(Hxx@X-Hx)/max(np.linalg.norm(Hx),1e-30)))
        else:X=None;S=None
        blocks.append(dict(H=H,Hxx=Hxx,Hxs=Hx,Hss=Hss,X=X,S=S,mechanical_mode=ux[:,0],threshold=threshold));eigenvalues.append(record)
    np.savez_compressed(folder/'CONDENSED_PRIMARY.npz',H=blocks[1]['H'],Hxx=blocks[1]['Hxx'],Hxs=blocks[1]['Hxs'],Hss=blocks[1]['Hss'],
        S=blocks[1]['S'] if blocks[1]['S'] is not None else np.empty((0,0)),X=blocks[1]['X'] if blocks[1]['X'] is not None else np.empty((p.nx,0)),critical=critical)
    save(folder/'FACTORIZATION.json',eigenvalues)
    lambdas=np.array([r['lambda_x_min'] for r in eigenvalues]);hxx_convergence=float(np.ptp(lambdas)/max(abs(lambdas[1]),1e-30))
    matrix_pass=all(r['Hxx_raw_asymmetry']<=cfg['raw_symmetry_tolerance'] and r['cross_relative_error']<=cfg['cross_symmetry_tolerance'] and r['Hss_relative_error']<=cfg['hss_tolerance'] for r in assembly)
    hss_pass=all(r['Hss_relative_error']<=cfg['hss_tolerance'] for r in assembly)
    if not hss_pass:raise RuntimeError('HSS_ANALYTIC_GATE_FAILED')
    oracle=mechanical_oracle(p,z,blocks[1]['Hxx']);save(folder/'MECHANICAL_ORACLE.json',oracle)
    hxx_pass=spec['smooth_restricted_energy'] and matrix_pass and oracle['passed'] and hxx_convergence<=cfg['FD_magnitude_relative_tolerance'] and len(set(np.sign(lambdas)))==1
    status='STABILITY_UNRESOLVED';mechanical='UNRESOLVED';selection=None;mode=None;direct=None;rayleigh=None;schur_tests=[]
    if hxx_pass:
        if all(r['lambda_x_min']<-r['curvature_tol'] for r in eigenvalues):
            mechanical='UNSTABLE';mode=np.r_[blocks[1]['mechanical_mode'],np.zeros(p.ns)];selection=dict(method='UNSTABLE_MECHANICAL_MODE');rayleigh=eigenvalues[1]['lambda_x_min']
        elif all(r['lambda_x_min']>r['curvature_tol'] for r in eigenvalues):
            mechanical='STABLE'
            if m:
                (rayleigh,name,pvec),selection=select_cone(blocks[1]['S'],groups,critical,s[critical]-p.lower[critical]);selection['selected_name']=name
                u=-blocks[1]['X']@pvec;mode=np.zeros_like(z);mode[:p.nx]=u;mode[p.nx+critical]=pvec;mode/=np.linalg.norm(mode)
                # Algebraic and physical-energy checks of the actual condensed
                # displacement for deterministic phase directions.
                rng=np.random.default_rng(cfg['random_seed']+1)
                for i in range(3):
                    vv=rng.normal(size=m)
                    vv[groups['WEAK_LOWER'][critical]]=abs(vv[groups['WEAK_LOWER'][critical]])
                    vv[groups['WEAK_UPPER'][critical]]=-abs(vv[groups['WEAK_UPPER'][critical]])
                    vv/=np.linalg.norm(vv);uu=-blocks[1]['X']@vv;w=np.r_[uu,vv];lhs=float(w@blocks[1]['H']@w);rhs=float(vv@blocks[1]['S']@vv)
                    full=np.zeros_like(z);full[:p.nx]=uu;full[p.nx+critical]=vv;full/=np.linalg.norm(full)
                    curv=lhs/(w@w);probe=direct_probe(p,z,full,curv)
                    schur_tests.append(dict(direction=i,algebraic_relative_error=abs(lhs-rhs)/max(abs(rhs),1e-30),direct=probe))
            else:mode=np.r_[blocks[1]['mechanical_mode'],np.zeros(p.ns)];selection=dict(method='STRONG_BOUND_MECHANICAL_CERTIFICATE');rayleigh=eigenvalues[1]['lambda_x_min']
        else:mechanical='MARGINAL';status='STATIONARY_MARGINAL'
    if mode is not None:
        keep=np.r_[np.arange(p.nx),p.nx+critical];v=mode[keep];curvatures=[float(v@b['H']@v) for b in blocks]
        direct=direct_probe(p,z,mode,curvatures[1]);feasible=cone_feasible(mode[p.nx+critical],groups,critical) and np.all(mode[p.nx:][groups['STRONG_LOWER']|groups['STRONG_UPPER']]==0)
        if direct['passed'] and feasible:
            if all(c<-b['threshold'] for c,b in zip(curvatures,blocks)):status='STATIONARY_UNSTABLE'
            elif mechanical=='STABLE' and (not m or all(r['lambda_S_min']>r['curvature_tol'] for r in eigenvalues)):
                status='STATIONARY_STABLE'
        direct.update(curvatures_at_three_FD_steps=curvatures,cone_feasible=bool(feasible))
        np.savez_compressed(state_folder/'MODE.npz',vector=mode,critical=critical,x=x,s=s,lower=p.lower)
        save(state_folder/'MODE_ENERGY.json',direct)
    result=dict(energy_oracle_revision=2,label=label,strain=p.strain,state_hash=folder.name,cache=str(folder),metrics=p.metrics(x,s),counts={k:int(v.sum()) for k,v in groups.items()},m=m,
        critical_IDs=critical.tolist(),mechanical_status=mechanical,status=status,spectral=spec,assembly=assembly,eigenvalues=eigenvalues,
        Hxx_gate=bool(hxx_pass),matrix_gate=bool(matrix_pass),Hss_gate=bool(hss_pass),Hxx_FD_relative_spread=hxx_convergence,
        mechanical_direct_curvature_max_error=oracle['max_relative_error'],lambda_x_min=eigenvalues[1]['lambda_x_min'],lambda_S_min=eigenvalues[1].get('lambda_S_min'),
        best_feasible_Rayleigh=rayleigh,mode_selection=selection,mode_direct_pass=bool(direct and direct['passed']),schur_tests=schur_tests,
        Hessian_assembly_seconds=sum(r['wall_seconds'] for r in assembly),new_Hessian_assembly_seconds=sum(r['wall_seconds'] for r in assembly if not r['cache_hit']),wall_seconds=time.perf_counter()-start,
        parent_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,worker_peak_RSS_bytes=max(r['worker_peak_RSS_bytes'] for r in assembly))
    save(state_folder/'STABILITY.json',result);print('STABILITY',label,'m',m,'Hxx',result['lambda_x_min'],'S',result['lambda_S_min'],status,'oracle',oracle['max_relative_error'],flush=True)
    return result


def existing():
    assert read(OUT/'INPUT_VERIFICATION.json')['passed'];cfg=settings();rows=[]
    with ProcessPoolExecutor(max_workers=cfg['workers'],mp_context=mp.get_context('spawn')) as pool:
        for label,e,source,summary in saved_inputs():
            data=old.ca.read_npz(OUT/f'states/{label}/STATE.npz');p=problem(e,data['lower']);key=state_key(p,data['x'],data['s']);path=OUT/f'states/{label}/STABILITY.json'
            if path.exists() and read(path)['state_hash']==key and read(path).get('energy_oracle_revision')==2:row=read(path)
            else:row=audit_state(label,p,data['x'],data['s'],pool)
            rows.append(row);save(OUT/'EXISTING_STATES.json',dict(states=rows,complete=len(rows)==9))


def new_stationary(label,strain,lower_row,pool):
    data=old.ca.read_npz(OUT/f'states/{lower_row["label"]}/STATE.npz');pstat=old.Problem(5,strain,data['s']);x0=old.tracker.displacement_predictor(pstat,data['x'],lower_row['strain']);z0=pstat.pack(x0,data['s'])
    folder=OUT/f'continuation/{label}'
    if (folder/'SUMMARY.json').exists():
        result=read(folder/'SUMMARY.json');a=old.ca.read_npz(folder/'STATE.npz')
        assert result['strain']==strain and np.array_equal(a['lower'],data['s']) and np.array_equal(a['z_initial'],z0)
        z=a['z'];print('REUSE STATIONARY',label,flush=True)
    else:result,z=petsc.solve(pstat,z0,folder,pool=pool,workers=settings()['workers'])
    result.update(label=label,previous_label=lower_row['label'],previous_strain=lower_row['strain']);save(folder/'SUMMARY.json',result)
    if not result['scientifically_accepted']:return None,result
    x,s=pstat.unpack(z);statefolder=OUT/f'states/{label}';statefolder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(statefolder/'STATE.npz',x=x,s=s,lower=pstat.lower,strain=strain)
    p=problem(strain,pstat.lower);row=audit_state(label,p,x,s,pool)
    return row,result


def localization(row):
    d=old.ca.read_npz(OUT/f'states/{row["label"]}/STATE.npz');v=old.ca.read_npz(OUT/f'states/{row["label"]}/MODE.npz')['vector'];p=problem(row['strain'],d['lower']);phase=abs(v[p.nx:]);mid=p.family.X[p.family.pairs].mean(axis=1)
    u=np.zeros_like(p.family.X);u[p.free]=v[:p.nx].reshape(-1,3)*p.length_scale;particle=old.v2.aggregate(p.family,phase)
    metrics=old.v2.participation(phase);center=np.average(mid,weights=phase,axis=0) if phase.sum() else np.zeros(3)
    metrics.update(center_m=center.tolist(),center_distance_to_clamp_m=float(old.v2.L/2-abs(center[2])),
        amplitude_weighted_clamp_distance_m=float(np.average(old.v2.L/2-abs(mid[:,2]),weights=phase)) if phase.sum() else None,
        top10_definition='Largest ceil(0.1*N) bond amplitudes divided by all bond amplitude sum',kind='NUMERICAL_VARIATIONAL_MODE_NOT_PHYSICAL_CRACK')
    np.savez_compressed(OUT/'CRITICAL_MODE.npz',X=p.family.X,midpoints=mid,phase=phase,particle_phase=particle,mechanical_mode=u)
    save(OUT/'MODE_LOCALIZATION.json',metrics);return metrics


def branch():
    existing_result=read(OUT/'EXISTING_STATES.json');assert existing_result['complete'];rows=list(existing_result['states']);cfg=settings();c=cfg['continuation'];bracket=None;reason=None;history=[]
    unstable=[r for r in rows if r['status']=='STATIONARY_UNSTABLE'];stable=[r for r in rows if r['status']=='STATIONARY_STABLE']
    pool=ProcessPoolExecutor(max_workers=cfg['workers'],mp_context=mp.get_context('spawn'))
    try:
        if unstable:
            upper=min(unstable,key=lambda r:r['strain']);lowers=[r for r in stable if r['strain']<upper['strain']];lower=max(lowers,key=lambda r:r['strain']) if lowers else None
            if lower:bracket=(lower,upper)
        elif all(r['status']=='STATIONARY_STABLE' for r in rows):
            previous=rows[-1];anchor=previous['strain'];k=0
            while True:
                if previous['metrics']['maximum_s']>=c['max_s']:reason='MAX_PHASE_AUDIT_WINDOW_REACHED';break
                k+=1;strain=anchor+k*c['step']
                if strain>c['max_strain']:reason='MAX_STRAIN_AUDIT_WINDOW_REACHED';break
                row,result=new_stationary(f'E{k:02d}',strain,previous,pool)
                if row is None:reason='NEW_STATIONARY_STATE_REJECTED';break
                rows.append(row);save(OUT/'ALL_STATES.json',dict(states=rows))
                if row['status']=='STATIONARY_UNSTABLE':bracket=(previous,row);reason='FIRST_VERIFIED_UNSTABLE';break
                if row['status']!='STATIONARY_STABLE':reason='STABILITY_UNRESOLVED_STOP';break
                previous=row
        else:reason='EXISTING_STABILITY_UNRESOLVED_NO_EXTENSION'
        if bracket:
            lower,upper=bracket
            for k in range(c['max_refinements']):
                if upper['strain']-lower['strain']<=c['refinement_width']:break
                mid=(lower['strain']+upper['strain'])/2;row,result=new_stationary(f'B{k:02d}',mid,lower,pool)
                if row is None:reason='BRACKET_STATIONARITY_UNRESOLVED';break
                rows.append(row);history.append(dict(lower_before=lower['label'],upper_before=upper['label'],midpoint=row['label'],strain=mid,status=row['status']))
                if row['status']=='STATIONARY_STABLE':lower=row
                elif row['status']=='STATIONARY_UNSTABLE':upper=row
                else:reason='NONSMOOTH_OR_NUMERICAL_UNRESOLVED';break
            nonsmooth=lower['counts']!=upper['counts'] or lower['spectral']['near_zero_BAQP']!=upper['spectral']['near_zero_BAQP']
            b=dict(interval=[lower['strain'],upper['strain']],width=upper['strain']-lower['strain'],lower_label=lower['label'],upper_label=upper['label'],history=history,
                kind='NONSMOOTH_STABILITY_BRACKET' if nonsmooth else 'STABILITY_BRACKET',resolved=upper['strain']-lower['strain']<=c['refinement_width'])
            save(OUT/'STABILITY_BRACKET.json',b);localization(upper)
        else:save(OUT/'STABILITY_BRACKET.json',dict(interval=None,width=None,history=[],resolved=False,reason=reason))
        save(OUT/'ALL_STATES.json',dict(states=rows,stop_reason=reason,coarse_step=c['step'],existing_roots_resolved_before_extension=True))
    finally:pool.shutdown()


if __name__=='__main__':
    import sys
    if sys.argv[1]=='prepare':declare();reproduce()
    elif sys.argv[1]=='existing':existing()
    elif sys.argv[1]=='branch':branch()

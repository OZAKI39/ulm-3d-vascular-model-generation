"""V2 orchestration: unchanged stationary solver, dense C3 reference, validated HVP."""
from scripts.spatial_audit_v2 import *
from scripts.run_monolithic_stability_audit import full_hessian
from scipy.linalg import eigh,solve
from scipy.sparse.linalg import eigsh,cg,minres
from scipy.interpolate import LinearNDInterpolator,NearestNDInterpolator,RegularGridInterpolator

def dense_reference(p,z,label):
    folder=OUT/f'matrix_free/{label}';folder.mkdir(parents=True,exist_ok=True);blocks=[];rows=[]
    x,s=p.unpack(z);expected=2*p.state(x,s).W0*p.weights/p.energy_scale
    for eps in EPS:
        H,schemes=full_hessian(p,z,eps);S=(H+H.T)/2;xx=S[:p.nx,:p.nx];xs=S[:p.nx,p.nx:]
        ev=eigh(xx,eigvals_only=True);schur=S[p.nx:,p.nx:]-xs.T@solve(xx,xs,assume_a='pos')
        row=dict(epsilon=eps,raw_symmetry_error=float(np.linalg.norm(H-H.T)/np.linalg.norm(H)),
            Hss_analytic_error=float(np.linalg.norm(H[p.nx:,p.nx:]-np.diag(expected))/np.linalg.norm(expected)),
            Hxx_condition=float(ev[-1]/ev[0]),Hxx_min=float(ev[0]),Hxx_max=float(ev[-1]),Schur_min=float(eigh(schur,eigvals_only=True)[0]),
            Schur_condensation_residual=float(np.linalg.norm(xx@solve(xx,xs)-xs)/np.linalg.norm(xs)),
            curvature_threshold=float(1e-7*max(np.linalg.norm(S,2),1)))
        np.savez_compressed(folder/f'DENSE_{eps:g}.npz',H_raw=H,H_symmetric=S,Hxx=xx,Schur=schur)
        rows.append(row);blocks.append(S)
    save(folder/'DENSE_REFERENCE.json',dict(label='CORRECTED_C3_DENSE_REFERENCE',length_scale_m=p.length_scale,energy_scale_J=p.energy_scale,rows=rows))
    return blocks,rows

def validate_matrix_free():
    p=Problem(3,.01,common_scale=True);x=affine(p);z=p.pack(x,np.zeros(p.ns));blocks,rows=dense_reference(p,z,'C3_validation')
    H=blocks[1];_,ev=eigh(H);rng=np.random.default_rng(84621)
    vectors=[(f'random_{i}',rng.normal(size=len(z))) for i in range(10)]
    vectors += [('lowest_coupled_mode',ev[:,0]),('uniform_phase_opening',np.r_[np.zeros(p.nx),np.ones(p.ns)]),
        ('random_phase_opening',np.r_[np.zeros(p.nx),abs(rng.normal(size=p.ns))])]
    tests=[]
    for eps,HH in zip(EPS,blocks):
        op=hessian_operator(p,z,eps)
        for name,v in vectors:
            v=v/np.linalg.norm(v);a=HH@v;b=op@v
            tests.append(dict(epsilon=eps,direction=name,relative_error=float(np.linalg.norm(a-b)/np.linalg.norm(a))))
    # The operator chain rule must also reproduce the verified tracker Jacobian.
    with patch.object(tracker,'settings',lambda:read(OUT/'provenance/PREDECLARED_CONFIG.json')):
        J=tracker.residual_jacobian(p,z);M=matrix_free_residual_jacobian(p,z)
        jt=[]
        for _,v in vectors[:10]:
            v=v/np.linalg.norm(v)
            jt.append(dict(forward=float(np.linalg.norm(J@v-M@v)/np.linalg.norm(J@v)),
                transpose=float(np.linalg.norm(J.T@v-M.rmatvec(v))/np.linalg.norm(J.T@v))))
    passed=max(r['relative_error'] for r in tests)<5e-5 and max(max(r.values()) for r in jt)<5e-5
    result=dict(passed=passed,tolerance=5e-5,predeclared=True,tests=tests,stationary_Jacobian_tests=jt,
        max_relative_error=max(r['relative_error'] for r in tests),dense_reference_valid=all(r['raw_symmetry_error']<1e-5 and r['Hss_analytic_error']<1e-6 for r in rows))
    save(OUT/'matrix_free/C3_DENSE_VS_HVP.json',result);print('HVP validation',result['max_relative_error'],passed,flush=True)
    assert passed and result['dense_reference_valid']

def phase_directions(p,H=None,op=None):
    f=p.family;mid=f.X[f.pairs].mean(axis=1);state=p.state(affine(p),p.lower);rng=np.random.default_rng(84621)
    ds=[('A_uniform',np.ones(p.ns)),('B_driver',np.maximum(state.W0,0))]
    if H is not None:_,v=eigh(H,subset_by_index=[0,0]);v=v[:,0]
    else:_,vv=eigsh(op,k=1,which='SA',tol=1e-7,v0=np.random.default_rng(193).normal(size=p.nx+p.ns));v=vv[:,0]
    ds.extend([('C_positive_lowest',np.maximum(v[p.nx:],0)),('D_sign_reversed_lowest',np.maximum(-v[p.nx:],0))])
    source=OUT/'modes/C3_SELECTED.npz'
    if source.exists():
        data=ca.read_npz(source);interpolator=NearestNDInterpolator(data['midpoints'],abs(data['phase']))
        ds.append(('E_C3_midpoint_projection',interpolator(mid)))
    else:ds.append(('E_C3_self_midpoint_projection',np.maximum(v[p.nx:],0)))
    for i in range(12):ds.append((f'F_nonnegative_random_{i}',rng.random(p.ns)))
    for i,z in enumerate(np.unique(np.round(mid[:,2],14))):ds.append((f'G_layer_{i}',np.isclose(mid[:,2],z,atol=1e-14,rtol=0).astype(float)))
    for label,mask in [('clamp',abs(mid[:,2])>=L/4),('lateral',np.max(abs(mid[:,:2]),axis=1)>=L/4),('interior',np.max(abs(mid),axis=1)<L/4)]:
        ds.append(('H_'+label,mask.astype(float)))
    return [(name,q/np.linalg.norm(q)) for name,q in ds if np.linalg.norm(q)>1e-20]

def cone_probes(p,z,blocks,epsrows):
    H=blocks[1];x,s=p.unpack(z);masks,active=tracker.active_sets(p,x,s);strong=masks['LOWER_STRONG_ACTIVE'];weak=masks['LOWER_WEAK_NEAR_SWITCH']
    grad=p.objective_gradient(z)[1];tests=[];vectors=[];thresholds=[r['curvature_threshold'] for r in epsrows]
    for name,q in phase_directions(p,H=H):
        q=q.copy();q[strong]=0
        if np.linalg.norm(q)<1e-20:continue
        u=-solve(H[:p.nx,:p.nx],H[:p.nx,p.nx:]@q,assume_a='pos');v=np.r_[u,q];v/=np.linalg.norm(v)
        probes=[tracker.probe_point(p,z,v,a) for a in [1e-8,1e-7,1e-6]]
        curvature=[float(v@S@v) for S in blocks];mf=[float(v@(hessian_operator(p,z,e)@v)) for e in EPS]
        linear=float(grad@v);critical=bool(np.all(abs(v[p.nx:][strong])<1e-14) and np.all(v[p.nx:][weak]>=-1e-14) and abs(linear)<1e-8)
        negative=critical and any(t['feasible'] for t in probes) and all(r<-t for r,t in zip(curvature,thresholds)) and all(r<-t for r,t in zip(mf,thresholds))
        tests.append(dict(index=len(vectors),category=name,linear_variation=linear,critical_cone=critical,
            phase_feasible=bool(np.all(v[p.nx:][weak]>=0)),probes=probes,curvatures=curvature,matrix_free_curvatures=mf,
            condensation_relative_residual=float(np.linalg.norm(H[:p.nx,:p.nx]@u+H[:p.nx,p.nx:]@q)/max(np.linalg.norm(H[:p.nx,p.nx:]@q),1e-30)),
            verified_negative=negative));vectors.append(v)
    return tests,vectors,active

def mode_data(p,v,label):
    f=p.family;a=aggregate(f,v[p.nx:]);u=np.zeros_like(f.X);u[p.free]=v[:p.nx].reshape(-1,3)*p.length_scale
    mag=np.linalg.norm(u,axis=1);a/=max(a.max(),1e-30);mag/=max(mag.max(),1e-30)
    gridaxis=np.linspace(-L/2,L/2,21);grid=np.stack(np.meshgrid(gridaxis,gridaxis,gridaxis,indexing='ij'),axis=-1)
    axis=np.unique(f.X[:,0]);grid_a=RegularGridInterpolator((axis,axis,axis),a.reshape(p.n,p.n,p.n))(grid)
    grid_u=RegularGridInterpolator((axis,axis,axis),mag.reshape(p.n,p.n,p.n))(grid)
    center=np.average(f.X,weights=a,axis=0);xx=f.X-center;spread=np.einsum('i,ij,ik->jk',a,xx,xx)/a.sum()
    metrics=dict(**participation(a),center_m=center.tolist(),spread_tensor_m2=spread.tolist(),
        common_grid_shape=[21,21,21],clamp_distance_mean_m=float(np.average(L/2-abs(f.X[:,2]),weights=a)),
        lateral_distance_mean_m=float(np.average(L/2-np.max(abs(f.X[:,:2]),axis=1),weights=a)))
    masks,_=tracker.active_sets(p,*p.unpack(p.pack(affine(p),p.lower)))
    mid=f.X[f.pairs].mean(axis=1);layer=[]
    for zz in np.unique(f.X[:,2]):
        sel=f.X[:,2]==zz;bs=abs(mid[:,2]-zz)<=p.cloud.spacing/2+1e-15
        layer.append(dict(z_m=float(zz),weak_bond_fraction=float(masks['LOWER_WEAK_NEAR_SWITCH'][bs].mean()),
            interior_bond_fraction=float(masks['INTERIOR'][bs].mean()),mode_amplitude=float(a[sel].mean())))
    save(OUT/f'localization/{label}.json',dict(**metrics,layers=layer))
    np.savez_compressed(OUT/f'modes/{label}.npz',phase=v[p.nx:],displacement=u,particle_phase=a,particle_displacement=mag,
        grid_axis=gridaxis,grid_phase=grid_a,grid_displacement=grid_u,midpoints=mid,X=f.X)
    return metrics

def audit_state(n,strain,lower,x_initial,s_initial,label):
    start=time.time();p=Problem(n,strain,lower);res=stationary(p,x_initial,s_initial)
    folder=OUT/f'branch_tracking/{label}';folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/'STATE.npz',**{k:res[k] for k in ['x','s','z','z_initial','lower']})
    record={k:v for k,v in res.items() if k not in ['x','s','z','z_initial','lower']}
    record.update(level=f'C{n}',strain=strain,label=label,wall_seconds=time.time()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    save(folder/'STATIONARY.json',record)
    masks,active=tracker.active_sets(p,res['x'],res['s']);active['fractions']={k:v/p.ns for k,v in active['counts'].items()}
    active['bond_midpoints_m']=p.family.X[p.family.pairs].mean(axis=1).tolist();save(OUT/f'active_sets/{label}.json',active)
    state=p.state(res['x'],res['s']);reaction=float(-state.force[p.top,2].sum())
    summary=dict(label=label,strain=strain,stationary=res['accepted'],stationary_reason=res['reason'],gates=res['gates'],
        metrics=res['metrics'],active_counts=active['counts'],active_fractions=active['fractions'],reaction_N=reaction,
        nominal_stress_Pa=reaction/L**2,elastic_density_Pa=state.energy/L**3,
        crack_density_Pa=res['metrics']['crack_J']/L**3,total_density_Pa=res['metrics']['total_J']/L**3,
        crack_energy_per_area_J_m2=res['metrics']['crack_J']/L**2,
        status='STATIONARY_UNRESOLVED',negative_witness=None,localization=None)
    if res['accepted']:
        pp=Problem(n,strain,lower,common_scale=True);zz=pp.pack(res['x'],res['s'])
        if n==3:
            blocks,epsrows=dense_reference(pp,zz,label);tests,vectors,_=cone_probes(pp,zz,blocks,epsrows)
            valid=all(r['raw_symmetry_error']<1e-5 and r['Hss_analytic_error']<1e-6 for r in epsrows)
            negatives=[r for r in tests if r['verified_negative']]
            if valid and negatives:
                selected=min(negatives,key=lambda r:r['curvatures'][1]);summary['status']='STATIONARY_UNSTABLE';summary['negative_witness']=selected
                v=vectors[selected['index']];summary['localization']=mode_data(pp,v,label)
                np.savez_compressed(OUT/f'modes/C{n}_SELECTED.npz',phase=v[pp.nx:],midpoints=pp.family.X[pp.family.pairs].mean(axis=1),vector=v)
            elif valid and active['counts']['LOWER_STRONG_ACTIVE']==pp.ns and all(r['Hxx_min']>r['curvature_threshold'] for r in epsrows):
                summary['status']='STATIONARY_STABLE'
            else:
                keep=np.r_[np.arange(pp.nx),pp.nx+np.flatnonzero(~masks['LOWER_STRONG_ACTIVE'])]
                if valid and all(eigh(H[np.ix_(keep,keep)],eigvals_only=True,subset_by_index=[0,0])[0]>r['curvature_threshold'] for H,r in zip(blocks,epsrows)):
                    summary['status']='STATIONARY_STABLE'
            save(OUT/f'cone_stability/{label}.json',dict(candidates=tests,classification=summary['status'],dense=epsrows))
        else:
            # Fine stability is eligible only following a stationary accepted state.
            summary.update(status='STATIONARY_UNRESOLVED',stability_pending='MATRIX_FREE_SPECTRAL_AUDIT_REQUIRED')
    save(folder/'SUMMARY.json',summary)
    print(label,summary['status'],res['reason'], 'force=',res['metrics']['force_residual'], 'green=',res['metrics']['minimum_Green_strain'],flush=True)
    return res,summary

def schedule(n,name,start_strain=.009):
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json');ns=len(fixture(n)[1].pairs);lower=np.zeros(ns);p0=Problem(n,start_strain)
    xx=affine(p0);ss=lower.copy();prev_strain=start_strain;rows=[];last_stable=None;upper=None
    for k,e in enumerate(cfg['schedules'][name]):
        p=Problem(n,e,lower);xx=tracker.displacement_predictor(p,xx,prev_strain)
        result,row=audit_state(n,e,lower,xx,ss,f'C{n}_{name}_{k:03d}');rows.append(row)
        if row['status']=='STATIONARY_UNRESOLVED':break
        if row['status']=='STATIONARY_UNSTABLE':upper=(result,row);break
        last_stable=(result,row);xx=result['x'];ss=result['s'];lower=ss.copy();prev_strain=e
    initial_bracket=[last_stable[1]['strain'],upper[1]['strain']] if last_stable and upper else None
    if initial_bracket:
        for k in range(12):
            if upper[1]['strain']-last_stable[1]['strain']<=1e-5:break
            lo=last_stable;mid=(lo[1]['strain']+upper[1]['strain'])/2;p=Problem(n,mid,lo[0]['s'])
            xp=tracker.displacement_predictor(p,lo[0]['x'],lo[1]['strain'])
            result,row=audit_state(n,mid,lo[0]['s'],xp,lo[0]['s'],f'C{n}_{name}_refine_{k:02d}');rows.append(row)
            if row['status']=='STATIONARY_STABLE':last_stable=(result,row)
            elif row['status']=='STATIONARY_UNSTABLE':upper=(result,row)
            else:break
    valid=last_stable is not None and upper is not None and not any(r['status']=='STATIONARY_UNRESOLVED' and last_stable[1]['strain']<r['strain']<upper[1]['strain'] for r in rows)
    bracket=[last_stable[1]['strain'],upper[1]['strain']] if valid else None
    answer=dict(level=f'C{n}',schedule=name,initial_bracket=initial_bracket,interval=bracket,width=bracket[1]-bracket[0] if bracket else None,
        last_stable=last_stable[1] if last_stable else None,first_unstable=upper[1] if upper else None,
        complete=bool(valid and bracket[1]-bracket[0]<=1e-5),states=rows,unresolved_barrier=next((r for r in rows if r['status']=='STATIONARY_UNRESOLVED'),None))
    save(OUT/f'critical_brackets/C{n}_{name}.json',answer)
    return answer

def main():
    validate_matrix_free()
    all_results=[]
    for n in LEVELS:
        for name in ['A','B']:
            all_results.append(schedule(n,name));save(OUT/'critical_brackets/ALL.json',all_results)
    # A separate lower-load domain diagnostic does not bridge a rejected branch trial.
    for n in [5,7]:
        if not next(r for r in all_results if r['level']==f'C{n}')['complete']:
            for e in [0.,.00025]:
                p=Problem(n,e);audit_state(n,e,p.lower,affine(p),p.lower,f'C{n}_low_start_{e:g}')

if __name__=='__main__':main()

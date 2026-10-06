"""Frozen integration audit; run prepare, then exactly five serial solves."""
import argparse,csv,gc,hashlib,json,multiprocessing as mp,threading,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import psutil
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.higher_order_rk.mechanics_candidate import Reference,GradientOperator,VariationalBAQP
from pd_clot.pfpd_candidate.higher_order_rk.variational_energy import bond_state
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import BlockBoundary
from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_backend import PD2Backend,PD2SolverCallback
from pd_clot.pfpd_candidate.nonlocal_boundary.pd2_bulk_reference import PD2Stencil,oracle_job,field_for
from pd_clot.pfpd_candidate.nonlocal_boundary.localized_conversion import convert_localized
from pd_clot.pfpd_candidate.nonlocal_boundary.diagnostics import load_checks,support_profile
from pd_clot.pfpd_candidate.nonlocal_boundary.field_bias import bias_samples,bias_summary
from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_adapter import solve

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];OLD=ROOT/'verification/higher_order_baqp_rk_oracle_v1'
CFG=json.loads((OUT/'PD2_BOUNDARY_CONFIGURATION.json').read_text());H=CFG['h_m'];L=CFG['L_m'];DELTA=CFG['delta_m']
MAT=Material('svk',1000.,2000./3);KINDS=CFG['solve_order'][:4]
def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x))
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2,default=plain,allow_nan=False)+'\n')
def load(name):return json.loads((OUT/name).read_text())
def csvsave(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def norm(x):return np.linalg.norm(x,axis=tuple(range(1,x.ndim)))
def require(ok,reason):
    if not ok:save('HARD_STOP.json',dict(reason=reason));raise RuntimeError(reason)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def integrity():
    old=load('PROTECTED_INPUT_HASHES.json');changed=[p for p,h in old.items() if sha(ROOT/p)!=h]
    return dict(protected_files=len(old),changed=changed,PASS=not changed)
def primary():
    with np.load(OLD/'REFERENCE_GEOMETRY.npz') as a:
        f=Family.from_arrays(a['X'],a['volume'],a['pairs'],DELTA,a['omega']);s,t=a['source'],a['target']
        require(np.array_equal(Reference.from_family(f).m,a['bond_weight']),'Accepted bond weights changed')
    with np.load(OLD/'PD2_REFERENCE_CACHE.npz') as a:g=a['g']
    b=PD2Backend(VariationalBAQP(Reference.from_family(f),GradientOperator('PD2',s,t,g),MAT),f)
    return BlockBoundary(f.X.min(0),f.X.max(0)),f,b
def scales(kind,field):
    if 'affine' in kind and not kind.startswith('nonaffine'):
        p=np.linalg.norm(field.response()['P']);return dict(F_ref_N=p*H**2,solver_force_floor_N=p*L**2,gradient_scale=.003,b_ref_N_m3=0.)
    s=field.scales(H);s['solver_force_floor_N']=s['F_ref_N'];return s
def analytic(field,X):
    return np.broadcast_to(field.F,(len(X),3,3)) if hasattr(field,'F') else field.response(X[:,2])['F']
def pool(specs):
    rows=[]
    with ProcessPoolExecutor(CFG['oracle_workers'],mp_context=mp.get_context('spawn')) as p:
        for r in p.map(oracle_job,specs):
            print('ORACLE',r['kind'],r.get('label'),r.get('level'),r['q_bulk_N_m3'],flush=True);rows.append(r)
    return rows

def prepare():
    g,f,b=primary();saved={};summaries={};support=[];comparison=[];work=[]
    # Actual reversed force graph, including endpoint volumes and weight support.
    closures=[]
    for n in [19,25,27]:
        s=PD2Stencil(n);closures.append(s.closure(s.id_at_offset([0,0,0])))
    save('PD2_FORCE_STENCIL_CLOSURE.json',dict(domains=closures,selected_oracle_N=27,
        minimal_odd_for_all_required_bulk_volumes=27,explanation='D3 maximum coordinate offset is 12h; a physical boundary at 12h has half-volume, so the first sufficient odd scaffold is N27 (boundary 13h).'))
    require(closures[-1]['all_required_nodes_complete'],'N27 PD2 closure insufficient')
    rows=pool([dict(kind=k,label='affine_first_gate',Z=0.) for k in KINDS[:2]])
    affine=[]
    for r in rows:
        err=np.linalg.norm(r['q_bulk_N_m3'])*H**3/scales(r['kind'],field_for(r['kind']))['F_ref_N']
        affine.append(dict(kind=r['kind'],bulk_force_scaled=err,PASS=err<=1e-10))
    save('PD2_AFFINE_BULK_GATE.json',affine);require(all(r['PASS'] for r in affine),'Affine PD2 bulk force is not zero')
    zs=np.unique(f.X[:,2]);iz=np.searchsorted(zs,f.X[:,2])
    rows+=pool([dict(kind=k,label='primary',level=j,Z=float(z)) for k in KINDS[2:] for j,z in enumerate(zs)])
    rows+=pool([dict(kind=k,label='shifted_N29',level=j,Z=float(zs[j]),n=29,offset=[1,-1,1],point=[3*H,-2*H,float(zs[j])]) for k in KINDS[2:] for j in [0,9,18]])
    save('PD2_BULK_ORACLE_RAW.json',rows)
    bulkrows=[];validation=[]
    for k in KINDS:
        field=field_for(k);sc=scales(k,field)
        rr=[r for r in rows if r['kind']==k and r['label'] in ('primary','affine_first_gate')]
        q=np.array([r['q_bulk_N_m3'] for r in rr])
        if k in KINDS[:2]:q=np.broadcast_to(q[0],(19,3)).copy();body=np.zeros_like(q)
        else:body=field.response(zs)['b']
        for r in rows:
            if r['kind']==k and r['label']=='shifted_N29':
                e=np.linalg.norm(np.array(r['q_bulk_N_m3'])-q[r['level']])*H**3/sc['F_ref_N']
                validation.append(dict(kind=k,level=r['level'],shifted_error_scaled=e,PASS=e<=1e-10))
        for j,z in enumerate(zs):
            bulkrows.append(dict(case=k,level=j,Z_m=z,**{f'q_{a}_N_m3':q[j,c] for c,a in enumerate('xyz')},**{f'b_{a}_N_m3':body[j,c] for c,a in enumerate('xyz')},**{f'e_model_{a}_N_m3':(q+body)[j,c] for c,a in enumerate('xyz')}))
        x=field.position(f.X);st=b.evaluate(x)
        arr=convert_localized(b,x,q[iz],body[iz],f.volume)
        closure=float(norm(arr['f_FD']+arr['f_BC_localized']+arr['f_body_local']-arr['f_model']).max()/sc['F_ref_N'])
        target_res=float(norm(st.force+arr['f_external_NL']).max()/sc['F_ref_N'])
        model=q+body
        summaries[k]=dict(scales=sc,nodewise_decomposition_error=closure,NL_target_residual_scaled=target_res,
            bulk_model_max_N_m3=float(norm(model).max()),bulk_model_RMS_N_m3=float(np.sqrt(np.mean(norm(model)**2))),
            bulk_model_relative_to_body=float(norm(model).max()/sc['b_ref_N_m3']) if sc['b_ref_N_m3'] else None,target_energy_J=st.energy)
        require(max(closure,target_res)<=1e-10,'Nodewise decomposition failed')
        # Accepted prior candidate state regression (strict same arrays where x matches).
        previous={'normal_affine':'AFFINE_AXIAL','shear_affine':'AFFINE_SHEAR','nonaffine_axial':'AXIAL','nonaffine_shear':'SHEAR'}[k]
        with np.load(ROOT/f'verification/higher_order_baqp_energy_force_candidate_v1/{previous}_STATE.npz') as a:
            summaries[k]['previous_state_x_identical']=bool(np.array_equal(x,a['x']))
            if np.array_equal(x,a['x']):
                summaries[k]['previous_force_identical']=bool(np.array_equal(st.force,a['force']))
                summaries[k]['previous_energy_identical']=bool(st.energy==a['energy_J'])
                require(summaries[k]['previous_force_identical'] and summaries[k]['previous_energy_identical'],'Accepted PD2 state changed')
        saved.update({k+'__'+a:v for a,v in arr.items()});saved[k+'__target']=x;saved[k+'__q_bulk']=q[iz];saved[k+'__b_local']=body[iz]
        if k in KINDS[:2]:
            checks=load_checks(g,f.X,arr['f_BC_localized'],field)
            summaries[k]['load_checks']=checks;summaries[k]['face_targets']=field.exact_face_integrals(g)
            work.append(dict(case=k,channel='combined',force_error=checks['resultant_error'],moment_error=checks['moment_error'],work_error=checks['virtual_work_error']))
            require(all(checks[a] for a in ['resultant_PASS','moment_PASS','virtual_work_PASS']),'Affine force/moment/work gate failed')
            if k=='shear_affine':
                P=field.response()['P'];channels={}
                for channel,Pc in [('normal',P*np.eye(3)),('tangential',P-P*np.eye(3))]:
                    force=-b.stress_force(Pc);channels[channel]=force
                    ck=load_checks(g,f.X,force,field,P_override=Pc)
                    summaries[k][channel+'_load_checks']=ck
                    work.append(dict(case=k,channel=channel,force_error=ck['resultant_error'],moment_error=ck['moment_error'],work_error=ck['virtual_work_error']))
                    require(all(ck[a] for a in ['resultant_PASS','moment_PASS','virtual_work_PASS']),'Shear stress channel consistency failed')
                summaries[k]['stress_channel_sum_error']=float(norm(channels['normal']+channels['tangential']+st.force).max()/sc['F_ref_N'])
                require(summaries[k]['stress_channel_sum_error']<=1e-10,'Stress channel sum failed')
        oldpath={'normal_affine':'nonlocal_boundary_formulation_v1/NORMAL_AFFINE.npz','shear_affine':'nonlocal_boundary_formulation_v1/SHEAR_AFFINE.npz','nonaffine_axial':'nonlocal_boundary_nonaffine_localization_v1/NODEWISE_RESIDUAL_DECOMPOSITION.npz','nonaffine_shear':'nonlocal_boundary_nonaffine_shear_localization_v1/NODEWISE_SHEAR_DECOMPOSITION.npz'}[k]
        with np.load(ROOT/'verification'/oldpath) as a:
            require(np.array_equal(a['X'],f.X),'Historical C1 node mismatch')
            oldload=a['load' if k in KINDS[:2] else 'f_BC_localized']
        saved[k+'__C1_boundary']=oldload
        for backend,loadforce in [('C1',oldload),('PD2',arr['f_BC_localized'])]:
            pro,qu=support_profile(g,f,loadforce,H)
            support.extend(dict(case=k,backend=backend,depth_um=r['depth_m']*1e6,**r) for r in pro)
            row=dict(case=k,backend=backend,L1_N=float(norm(loadforce).sum()),L2_N=float(np.linalg.norm(loadforce)),max_N=float(norm(loadforce).max()))
            for qname,qr in qu.items():
                row.update({qname+'_um':qr['depth_m']*1e6,qname+'_h':qr['depth_h'],qname+'_delta':qr['depth_over_delta']})
            row.update(C1_minus_PD2_L2_N=float(np.linalg.norm(oldload-arr['f_BC_localized'])),C1_minus_PD2_max_N=float(norm(oldload-arr['f_BC_localized']).max()))
            comparison.append(row)
        # Full finite-domain vs selective graph (nonaffine force must agree).
        if k in KINDS[2:]:
            ss=PD2Stencil(19)
            for off in [[0,0,0],[-5,0,0],[0,0,8]]:
                pid=ss.id_at_offset(off);err=float(np.linalg.norm(ss.force(pid,field)-st.force[pid])/sc['F_ref_N'])
                validation.append(dict(kind=k,offset=off,selective_full_error_scaled=err,PASS=err<=1e-10))
            del ss
        print('PREPARED',k,'decomposition',closure,'model',summaries[k]['bulk_model_max_N_m3'],flush=True)
    save('PD2_ORACLE_VALIDATION.json',validation);require(all(r['PASS'] for r in validation),'PD2 oracle validation failed')
    csvsave('PD2_BULK_ORACLE.csv',bulkrows);csvsave('PD2_BOUNDARY_SUPPORT.csv',support);csvsave('C1_PD2_BOUNDARY_COMPARISON.csv',comparison);csvsave('PD2_FORCE_MOMENT_WORK.csv',work)
    np.savez_compressed(OUT/'PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz',X=f.X,volume=f.volume,depth_over_delta=g.distance(f.X)/DELTA,**saved)
    save('PD2_TARGET_SUMMARIES.json',summaries)
    # Independent translated N29 oracle vs N27 lateral-depth observations.
    depth=[]
    for k in KINDS[2:]:
        field=field_for(k);sc=scales(k,field);s=PD2Stencil(27)
        q=np.array(next(r['q_bulk_N_m3'] for r in rows if r['kind']==k and r['label']=='shifted_N29' and r['level']==9))
        for off in [[-9,0,0],[-5,0,0],[-1,0,0],[0,0,0]]:
            pid=s.id_at_offset(off);fd=s.force(pid,field);bc=q*s.V[pid]-fd;c=s.closure(pid)
            depth.append(dict(case=k,offset=off,depth_h=int(np.minimum(s.ijk[pid],26-s.ijk[pid]).min()),force_FD_N=fd,f_BC_N=bc,boundary_force_scaled=float(np.linalg.norm(bc)/sc['F_ref_N']),fully_closed=c['all_required_nodes_complete'],closure=c))
        require(depth[-1]['fully_closed'] and depth[-1]['boundary_force_scaled']<=1e-10,'Closed PD2 center localization failed')
    save('PD2_CLOSED_CENTER_LOCALIZATION.json',depth)
    # Solver adapter derivative independently checked, without equilibrium runs.
    cb=PD2SolverCallback(b,np.zeros_like(f.X));x=field_for('nonaffine_shear').position(f.X);cb(x)
    q=(f.X-g.lower)/g.lengths;v=H*np.prod(np.sin(np.pi*q),axis=1)[:,None]*np.array([.7,-.4,1.])
    hfd=1e-5;jv=cb.tangent(v);fd=(b.force(x+hfd*v)-b.force(x-hfd*v))/(2*hfd)
    err=float(np.linalg.norm(jv-fd)/np.linalg.norm(fd));w=v[:,[2,0,1]]
    sym=float(abs(np.sum(w*jv)-np.sum(v*cb.tangent(w)))/max(abs(np.sum(w*jv)),1e-30))
    save('PD2_SOLVER_ADAPTER_CHECK.json',dict(central_difference_step=hfd,JVP_relative_error=err,reciprocity_error=sym,PASS=max(err,sym)<=1e-7))
    require(max(err,sym)<=1e-7,'PD2 solver callback derivative check failed')
    save('PREPARATION_COMPLETE.json',dict(PASS=True))

def quality(g,f,b,x,kind):
    field=field_for(kind);scale=scales(kind,field)['gradient_scale'];ref=b.candidate.reference
    F=b.node_F(x);Ft=bond_state(ref,x,F,slice(None),MAT)[0];mid=.5*(f.X[f.pairs[:,0]]+f.X[f.pairs[:,1]])
    exact=analytic(field,f.X);base=bias_samples(exact,exact,MAT);psi_scale=float(base['Psi_plus_analytic'].max())
    summary={};profiles=[];saved={}
    for name,X,A,weight in [('Fbar',f.X,F,f.volume),('Ftilde',mid,Ft,ref.m)]:
        d=bias_samples(A,analytic(field,X),MAT);summary[name]=bias_summary(d,weight,scale,psi_scale)
        bins=np.rint(2*g.distance(X)/H).astype(int)
        for k in np.unique(bins):
            mask=bins==k;row=bias_summary({a:v[mask] for a,v in d.items()},weight[mask],scale,psi_scale)
            profiles.append(dict(case=kind,sample=name,depth_h=k/2,depth_over_delta=k/2*H/DELTA,**row))
        saved.update({name+'__'+a:v for a,v in d.items()});saved[name+'__depth_over_delta']=g.distance(X)/DELTA
    return summary,profiles,saved,psi_scale

def recovery(kind):
    require(load('PREPARATION_COMPLETE.json')['PASS'],'Preparation incomplete')
    require(not (OUT/(kind.upper()+'_RECOVERED_STATE.npz')).exists(),'Refusing to rerun an existing equilibrium case')
    g,f,b=primary();base='nonaffine_shear' if kind=='classical_shear' else kind;field=field_for(kind);sc=scales(base,field)
    with np.load(OUT/'PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz') as a:
        target=a[base+'__target'];external=a[base+('__f_external_classical' if kind=='classical_shear' else '__f_external_NL')]
        model=a[base+'__f_model']
    if kind=='classical_shear':
        e=float(norm(b.force(target)+external-model).max()/sc['F_ref_N'])
        save('PD2_CLASSICAL_TARGET_RESIDUAL.json',dict(error_scaled=e,model_force_max_N=float(norm(model).max()),PASS=e<=1e-10))
        require(e<=1e-10,'Classical target residual does not equal model V')
    fixed=g.dirichlet_mask(f.X);q=(f.X-g.lower)/g.lengths
    initial=target+.02*.003*L*np.prod(np.sin(np.pi*q),axis=1)[:,None]*np.array([1.,-.4,.7]);initial[fixed]=target[fixed]
    np.savez_compressed(OUT/(kind.upper()+'_SOLVE_INPUT.npz'),initial=initial,external=external,fixed=fixed)
    print('START EQUILIBRIUM',kind,'free dofs',3*(~fixed).sum(),flush=True)
    x,info=solve(f,MAT,target,fixed,initial,external,H,sc['solver_force_floor_N'],CFG,mechanics=b)
    st=b.evaluate(x);res=st.force+external;scale=max(float(np.linalg.norm(external[~fixed])),sc['solver_force_floor_N'])
    residual=float(np.linalg.norm(res[~fixed])/scale)
    summary,profiles,samples,psi=quality(g,f,b,x,base)
    for row in profiles:row['case']=kind
    np.savez_compressed(OUT/(kind.upper()+'_FIELD_DIAGNOSTICS.npz'),**samples)
    csvsave(kind.upper()+'_FIELD_PROFILE.csv',profiles)
    old=load('PD2_TARGET_SUMMARIES.json')[base]
    result=dict(case=kind,solver=info,force_residual=residual,solver_force_scale_N=scale,fields=summary,Psi_ref_Pa=psi,
        displacement_shift_max_m=float(norm(x-target).max()),initial_perturbation_max_m=float(norm(initial-target).max()),
        energy_target_J=old['target_energy_J'],energy_solved_J=st.energy,external_work_target_J=float(np.sum(external*(target-f.X))),
        external_work_solved_J=float(np.sum(external*(x-f.X))),work_interpretation='Work of fixed nodal load at this displacement; excludes grip loading-path work, not a conservative-path energy equality',
        fixed_displacement_error_m=float(norm(x[fixed]-target[fixed]).max()),
        PASS=bool(info['converged'] and residual<=1e-8 and (kind=='classical_shear' or all(v['F_error_max_scaled']<=1e-5 and v['Psi_plus_error_max_scaled']<=1e-4 for v in summary.values()))))
    save('PD2_'+kind.upper()+'.json',result)
    np.savez_compressed(OUT/(kind.upper()+'_RECOVERED_STATE.npz'),X=f.X,x=x,target=target,Fbar=st.Fbar,force=st.force,external=external,reaction=-res*fixed[:,None],fixed=fixed)
    print('DONE',kind,'residual',residual,'F errors',[v['F_error_max_scaled'] for v in summary.values()], 'PASS',result['PASS'],flush=True)
    require(result['PASS'],'Equilibrium or recovery gate failed: '+kind)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare']+CFG['solve_order']);a=ap.parse_args()
    require(integrity()['PASS'],'Protected files changed');require(psutil.virtual_memory().available/psutil.virtual_memory().total>=.2,'RAM headroom')
    start=time.monotonic();stop=threading.Event();p=psutil.Process();stats=dict(stage=a.stage,oracle_workers=8 if a.stage=='prepare' else 0,equilibrium_workers=0 if a.stage=='prepare' else 1,BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',peak_process_tree_RSS_bytes=0,minimum_available_RAM_fraction=1.)
    def monitor():
        while not stop.is_set():
            rss=0
            for proc in [p]+p.children(recursive=True):
                try:rss+=proc.memory_info().rss
                except psutil.Error:pass
            stats['peak_process_tree_RSS_bytes']=max(stats['peak_process_tree_RSS_bytes'],rss)
            stats['minimum_available_RAM_fraction']=min(stats['minimum_available_RAM_fraction'],psutil.virtual_memory().available/psutil.virtual_memory().total)
            stop.wait(.25)
    t=threading.Thread(target=monitor,daemon=True);t.start()
    try:
        prepare() if a.stage=='prepare' else recovery(a.stage)
    finally:
        stop.set();t.join();stats['wall_seconds']=time.monotonic()-start;save('RESOURCES_'+a.stage+'.json',stats);save('HISTORICAL_INTEGRITY.json',integrity())
    require(stats['minimum_available_RAM_fraction']>=.2,'RAM guard during stage')

if __name__=='__main__':main()

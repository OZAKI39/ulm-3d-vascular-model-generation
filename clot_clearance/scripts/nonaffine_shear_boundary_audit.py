"""Shear transferability audit using the existing nonaffine boundary framework."""
import argparse,csv,gc,hashlib,json,multiprocessing as mp,resource,threading,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import psutil
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.energy import Material,constitutive
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import block_family
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticShear,stress_channel
from pd_clot.pfpd_candidate.nonlocal_boundary.bulk_reference import ManufacturedStencil,oracle_job
from pd_clot.pfpd_candidate.nonlocal_boundary.localized_conversion import localize
from pd_clot.pfpd_candidate.nonlocal_boundary.residual_decomposition import balances
from pd_clot.pfpd_candidate.nonlocal_boundary.volume_force_conversion import stress_channel_load
from pd_clot.pfpd_candidate.nonlocal_boundary.diagnostics import support_profile
from pd_clot.pfpd_candidate.nonlocal_boundary.field_bias import bias_samples,bias_summary
from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_adapter import solve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/nonlocal_boundary_nonaffine_shear_localization_v1'
OLD=ROOT/'verification/nonlocal_boundary_formulation_v1'


def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x))
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,default=plain,allow_nan=False)+'\n')
def load(n):return json.loads((OUT/n).read_text())
def csvsave(n,rows):
    with (OUT/n).open('w',newline='') as s:
        w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def material():return Material(model='svk',mu=1000.,bulk=666.6666666666666)
def normmax(x):return float(np.linalg.norm(x,axis=-1).max())
def integrity():
    hashes=load('PROTECTED_INPUT_HASHES.json');changed=[p for p,h in hashes.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    return dict(protected_files=len(hashes),changed=changed,PASS=not changed)
def require(ok,why):
    if not ok:
        save('HARD_STOP.json',dict(reason=why,thresholds_not_changed=True));raise RuntimeError('HARD STOP: '+why)
def pool_jobs(specs,cfg):
    result=[]
    with ProcessPoolExecutor(cfg['oracle_workers'],mp_context=mp.get_context('spawn')) as pool:
        for row in pool.map(oracle_job,specs):
            result.append(row);print('ORACLE',row.get('label'),row.get('channel','combined'),row.get('level'),row['q_bulk_N_m3'],flush=True)
    return result


def reference_affine(cfg):
    field=QuadraticShear();h=cfg['h_m'];sc=field.scales(h);zs=np.linspace(-field.length/2,field.length/2,19)
    response=field.response(zs);W,P=constitutive(response['F'],material())
    constitutive_error=float(np.linalg.norm(P-response['P'],axis=(1,2)).max()/sc['P_ref_Pa'])
    energy_error=float(abs(W-response['W']).max()/max(response['W']))
    continuum={c:field.continuum(c) for c in ['combined','normal','tangential']}
    maximum=max(continuum['combined'][k] for k in ['force_closure_scaled','moment_closure_scaled','max_work_closure_scaled','current_moment_closure_scaled'])
    reference=dict(Z_m=zs,response=response,scales=sc,continuum=continuum,independent_constitutive_error=constitutive_error,
                   independent_energy_error=energy_error,VERIFIED=max(maximum,constitutive_error,energy_error)<=cfg['gate'],
                   formulas=dict(u_x='gamma Z + beta_s Z^2/L',k='gamma + 2 beta_s Z/L',P='mu [[k^2,0,k+k^3],[0,0,0],[k,0,k^2]]',
                                 W='mu*k^2/2+mu*k^4/4',b='[-mu(1+3k^2)k_prime,0,-2mu*k*k_prime]',traction='P N on all six faces',
                                 reference_moment='surface(X)+body(X)=integral sum_a e_a cross P_col_a, not generally zero'))
    save('EXACT_NONAFFINE_SHEAR_REFERENCE.json',reference);save('NONAFFINE_SHEAR_EXACT_REFERENCE.json',reference)
    require(reference['VERIFIED'],'independent shear continuum reference')
    closure=[]
    for n,offset in [(27,[0,0,0]),(29,[1,-1,1])]:
        s=ManufacturedStencil(field,n);c=s.closure(s.id_at_offset(offset));c.update(n=n,offset=offset);closure.append(c)
    save('BULK_STENCIL_CLOSURE.json',closure);require(all(c['FULL_FORCE_STENCIL_CLOSED'] for c in closure),'closed oracle scaffold')
    del s;gc.collect()
    rows=pool_jobs([dict(kind='shear',label='affine',beta=0.,level=k,Z=float(zs[k])) for k in [0,9,18]],cfg)
    save('AFFINE_SHEAR_ORACLE.json',rows)
    qerr=max(np.linalg.norm(r['q_bulk_N_m3'])*h**3/sc['F_ref_N'] for r in rows)
    require(qerr<=cfg['gate'],'affine shear bulk force regression')
    g,f=block_family(19,field.length,cfg['delta_m']);aff=QuadraticShear(beta_s=0.);st=evaluate(f,aff.position(f.X),material())
    q=np.asarray(rows[1]['q_bulk_N_m3']);new=q*f.volume[:,None]-st.force
    with np.load(OLD/'SHEAR_AFFINE.npz') as old:
        require(np.array_equal(old['X'],f.X),'same affine physical grid')
        err=float(np.linalg.norm(new-old['load'])/np.linalg.norm(old['load']))
    save('AFFINE_SHEAR_LIMIT_RECOVERY.json',dict(relative_L2_error=err,bulk_force_max_scaled=qerr,PASS=err<=cfg['gate']))
    require(err<=cfg['gate'],'affine shear B1 regression')
    print('REFERENCE/AFFINE VERIFIED',maximum,err,flush=True)


def oracle(cfg):
    require(load('AFFINE_SHEAR_LIMIT_RECOVERY.json')['PASS'],'affine shear first')
    field=QuadraticShear();zs=np.linspace(-field.length/2,field.length/2,19);sc=field.scales(cfg['h_m']);specs=[]
    for channel in ['combined','normal','tangential']:
        specs += [dict(kind='shear',label='primary',channel=channel,level=k,Z=float(z)) for k,z in enumerate(zs)]
    for k in [0,9,18]:
        specs += [dict(kind='shear',label='shifted',level=k,Z=float(zs[k]),n=29,offset=[1,-1,1],x=3*cfg['h_m'],y=-2*cfg['h_m']),
                  dict(kind='shear',label='reflected_x',level=k,Z=float(zs[k]),gamma=-.003,beta=-.001)]
    rows=pool_jobs(specs,cfg);save('BULK_SHEAR_ORACLE_RAW.json',rows)
    q={c:np.array([r['q_bulk_N_m3'] for r in rows if r.get('channel')==c]) for c in ['combined','normal','tangential']}
    checks=[]
    for r in rows[57:]:
        target=q['combined'][r['level']].copy()
        if r['label']=='reflected_x':target[0]*=-1
        err=float(np.linalg.norm(np.asarray(r['q_bulk_N_m3'])-target)*cfg['h_m']**3/sc['F_ref_N'])
        checks.append(dict(label=r['label'],level=r['level'],force_scaled_error=err,PASS=err<=cfg['gate']))
    cherr=normmax(q['normal']+q['tangential']-q['combined'])*cfg['h_m']**3/sc['F_ref_N']
    b=field.response(zs)['b'];e=q['combined']+b;trans=float(abs(q['combined'][:,1]).max()*cfg['h_m']**3/sc['F_ref_N'])
    comp_b=abs(b).max(axis=0);comp_e=abs(e).max(axis=0);ratios=[float(comp_e[j]/comp_b[j]) if comp_b[j]>1e-12*sc['b_ref_N_m3'] else None for j in range(3)]
    result=dict(checks=checks,channel_sum_error_scaled=cherr,y_force_scaled=trans,VERIFIED=all(c['PASS'] for c in checks) and max(cherr,trans)<=cfg['gate'],
                max_q_bulk_N_m3=normmax(q['combined']),max_b_local_N_m3=normmax(b),max_e_model_N_m3=normmax(e),
                e_model_component_max_N_m3=comp_e,b_local_component_max_N_m3=comp_b,component_relative_errors=ratios,
                relative_model_difference=normmax(e)/sc['b_ref_N_m3'])
    result['SHEAR_BULK_MODEL_DIFFERENCE']=('SMALL' if result['relative_model_difference']<=cfg['bulk_model_small_threshold'] else 'MATERIAL') if result['VERIFIED'] else 'INCONCLUSIVE'
    save('BULK_SHEAR_ORACLE_VALIDATION.json',result)
    csvsave('BULK_BAQP_SHEAR_REFERENCE.csv',[dict(channel=c,level=k,Z_m=z,**{f'q_{a}_N_m3':q[c][k,j] for j,a in enumerate('xyz')},
                    **{f'b_{a}_N_m3':field.body(zs,c)[k,j] for j,a in enumerate('xyz')},
                    **{f'e_{a}_N_m3':(q[c]+field.body(zs,c))[k,j] for j,a in enumerate('xyz')},
                    q_norm_N_m3=float(np.linalg.norm(q[c][k])),b_norm_N_m3=float(np.linalg.norm(field.body(zs,c)[k])),
                    e_norm_N_m3=float(np.linalg.norm((q[c]+field.body(zs,c))[k]))) for c in q for k,z in enumerate(zs)])
    require(result['VERIFIED'],'bulk shear oracle validation');print('BULK VERIFIED',result['max_e_model_N_m3'],flush=True)


def field_bias(g,f,st,field,cfg):
    i,j=f.pairs.T;mid=.5*(f.X[i]+f.X[j]);exact_node=field.response(f.X[:,2])['F']
    psi_scale=float(bias_samples(exact_node,exact_node,material())['Psi_plus_analytic'].max())
    summaries={};profiles=[];saved={};deterministic={}
    for name,X,F,weights in [('Fbar',f.X,st.Fbar,f.volume),('Ftilde',mid,st.Ftilde,2*f.wbar*f.volume[i]*f.volume[j])]:
        analytic=field.response(X[:,2])['F'];d=bias_samples(F,analytic,material());repeat=bias_samples(F,analytic,material())
        deterministic[name]=all(np.array_equal(d[k],repeat[k]) for k in d);del repeat
        summaries[name]=bias_summary(d,weights,cfg['gradient_scale'],psi_scale)
        depth=g.distance(X);bins=np.rint(2*depth/cfg['h_m']).astype(int)
        for k in np.unique(bins):
            mask=bins==k;metrics=bias_summary({a:v[mask] for a,v in d.items()},weights[mask],cfg['gradient_scale'],psi_scale)
            profiles.append(dict(sample=name,depth_h=k/2,depth_over_delta=k/2*cfg['h_m']/cfg['delta_m'],samples=int(mask.sum()),**metrics))
        saved.update({name+'_'+k:v for k,v in d.items()});saved[name+'_depth_over_delta']=depth/f.horizon
        saved[name+'_weight']=weights
        summaries[name]['discrete_gradient_sha256']=hashlib.sha256(F.tobytes()).hexdigest()
        del analytic,d;gc.collect()
    maximum=max(v['F_error_max_scaled'] for v in summaries.values())
    classification='LOW' if maximum<cfg['kinematic_low_threshold'] else ('MODERATE' if maximum<=cfg['kinematic_material_threshold'] else 'MATERIAL')
    save('SHEAR_FIELD_BIAS_SUMMARY.json',dict(samples=summaries,Psi_plus_reference_scale_Pa=psi_scale,deterministic=deterministic,
                                           KINEMATIC_BOUNDARY_BIAS=classification,PSI_PLUS_BOUNDARY_BIAS='REPORTED_NOT_ACCEPTANCE_RATED',
                                           Ftilde_reference='Exact analytic F at reference bond midpoint; depth evaluated there',
                                           field_bias_is_diagnostic_only=True))
    csvsave('SHEAR_FIELD_BIAS.csv',profiles);np.savez_compressed(OUT/'SHEAR_FIELD_BIAS_SAMPLES.npz',**saved)
    return summaries


def decompose(cfg):
    require(load('BULK_SHEAR_ORACLE_VALIDATION.json')['VERIFIED'],'bulk oracle first')
    field=QuadraticShear();sc=field.scales(cfg['h_m']);Fref=sc['F_ref_N'];g,f=block_family(19,field.length,cfg['delta_m'])
    target=field.position(f.X);st=evaluate(f,target,material());zs,iz=np.unique(f.X[:,2],return_inverse=True)
    raw=load('BULK_SHEAR_ORACLE_RAW.json');q={c:np.array([r['q_bulk_N_m3'] for r in raw if r.get('channel')==c])[iz] for c in ['combined','normal','tangential']}
    nodal={};channel_results={};workrows=[];force_channels={}
    for c in ['combined','normal','tangential']:
        fd=st.force if c=='combined' else -stress_channel_load(f,st,stress_channel(st.P,c))
        a=localize(fd,q[c],field.body(f.X[:,2],c),f.volume);nodal[c]=a;force_channels[c]=fd
        continuum=field.continuum(c)
        force,moment,work=balances(f.X,a,continuum,Fref,field.length,
                                  variations=[(n,u) for n,u,_ in field.variations(f.X)],moment_source=continuum['reference_moment_source_Nm'])
        channel_results[c]=dict(force=force,moment=moment)
        workrows += [dict(channel=c,**r) for r in work]
    a=nodal['combined'];b=field.body(f.X[:,2]);closure=normmax(a['f_FD']+a['f_BC_localized']+a['f_body_local']-a['f_model'])/Fref
    naive=normmax(a['f_BC_naive']-a['f_BC_localized']+a['f_bulk_target'])/Fref;nl=normmax(st.force+a['f_external_NL'])/Fref
    cherr={k:normmax(nodal['normal'][k]+nodal['tangential'][k]-a[k])/Fref for k in a}
    selective=ManufacturedStencil(field,19);checks=[]
    for off in [[0,0,0],[-5,0,0],[0,0,8]]:
        pid=selective.id_at_offset(off);checks.append(dict(offset=off,error_scaled=float(np.linalg.norm(selective.force(pid)-st.force[pid])/Fref)))
    save('SELECTIVE_VS_FULL_SHEAR.json',checks);require(max(v['error_scaled'] for v in checks)<=cfg['gate'],'nonaffine shear selective/full')
    del selective;gc.collect()
    save('SHEAR_FORCE_BALANCE.json',channel_results['combined']['force'])
    moment=channel_results['combined']['moment'];moment['current_position_internal_moment_Nm']=np.cross(target,st.force).sum(axis=0)
    moment['current_internal_moment_scaled']=float(np.linalg.norm(moment['current_position_internal_moment_Nm'])/(Fref*field.length))
    save('SHEAR_MOMENT_BALANCE.json',moment);save('SHEAR_CHANNEL_ACCOUNTING.json',dict(channels=channel_results,nodewise_channel_sum_errors=cherr,
               interpretation='Diagonal/off-diagonal of actual discrete Piola stress; no unique physical-face ownership of distributed nodal loads'))
    csvsave('SHEAR_VIRTUAL_WORK.csv',workrows)
    identities=dict(closure_scaled=closure,naive_scaled=naive,NL_target_residual_scaled=nl,max_channel_sum_scaled=max(cherr.values()),
                    max_work_closure_scaled=max(v['closure_scaled'] for v in workrows),target_energy_J=st.energy,
                    PASS=max(closure,naive,nl,max(cherr.values()),max(v['closure_scaled'] for v in workrows))<=cfg['gate'])
    save('SHEAR_NODEWISE_IDENTITIES.json',identities);require(identities['PASS'],'shear nodewise/channel/work closure')
    np.savez_compressed(OUT/'NODEWISE_SHEAR_DECOMPOSITION.npz',X=f.X,volume=f.volume,x_target=target,Fbar_target=st.Fbar,
                        F_analytic=field.response(f.X[:,2])['F'],q_bulk=q['combined'],b_local=b,e_model=q['combined']+b,
                        distance_faces=g.face_distances(f.X),d_min_over_delta=g.distance(f.X)/f.horizon,d_min_over_h=g.distance(f.X)/cfg['h_m'],**a,
                        f_BC_normal=nodal['normal']['f_BC_localized'],f_BC_tangential=nodal['tangential']['f_BC_localized'])
    profiles=[];quantiles={}
    for name in ['localized','naive']:
        rows,qs=support_profile(g,f,a['f_BC_'+name],cfg['h_m'])
        for v in qs.values():v['depth_um']=v['depth_m']*1e6
        quantiles[name]=qs;profiles += [dict(correction=name,**v) for v in rows]
    csvsave('SHEAR_BOUNDARY_LOCALIZATION.csv',profiles);save('SHEAR_SUPPORT_QUANTILES.json',quantiles)
    bias=field_bias(g,f,st,field,cfg)
    print('DECOMPOSITION VERIFIED',closure,'BIAS',{k:v['F_error_max_scaled'] for k,v in bias.items()},flush=True)


def n27(cfg):
    require(load('SHEAR_NODEWISE_IDENTITIES.json')['PASS'],'nodewise first')
    field=QuadraticShear();s=ManufacturedStencil(field,27);Fref=field.scales(cfg['h_m'])['F_ref_N']
    reference=oracle_job(dict(kind='shear',label='independent_N29',Z=0.,n=29,offset=[1,-1,1]));q=np.asarray(reference['q_bulk_N_m3']);rows=[];closures=[]
    for off in cfg['N27_offsets']:
        pid=s.id_at_offset(off);fd=s.force(pid);bc=q*s.V[pid]-fd;c=s.closure(pid);closures.append(c)
        depth=int(np.minimum(s.ijk[pid],26-s.ijk[pid]).min())
        rows.append(dict(offset=off,point_X_m=s.X[pid],depth_h=depth,depth_over_delta=depth*s.h/s.delta,force_FD_N=fd,
                         bulk_target_N=q*s.V[pid],localized_N=bc,localized_norm_N=float(np.linalg.norm(bc)),localized_scaled=float(np.linalg.norm(bc)/Fref),
                         naive_norm_N=float(np.linalg.norm(fd)),naive_scaled=float(np.linalg.norm(fd)/Fref),fully_closed=c['FULL_FORCE_STENCIL_CLOSED']))
    result=dict(points=rows,independent_reference=reference,closure=closures,PASS=rows[-1]['fully_closed'] and rows[-1]['localized_scaled']<=cfg['gate'])
    save('N27_SHEAR_DEPTH_CHECK.json',result);require(result['PASS'],'N27 closed-center shear correction')
    print('N27 VERIFIED',[(v['depth_h'],v['localized_scaled']) for v in rows],flush=True)


def recovery(cfg,kind):
    require(load('N27_SHEAR_DEPTH_CHECK.json')['PASS'],'N27 first')
    if kind=='classical':require(load('SHEAR_NL_MANUFACTURED_RECOVERY.json')['PASS'],'NL recovery first')
    field=QuadraticShear();h=cfg['h_m'];sc=field.scales(h);g,f=block_family(19,field.length,cfg['delta_m'])
    with np.load(OUT/'NODEWISE_SHEAR_DECOMPOSITION.npz') as a:
        target=a['x_target'].copy();Ftarget=a['Fbar_target'].copy();Fan=a['F_analytic'].copy()
        external=a['f_external_NL' if kind=='NL' else 'f_external_classical'].copy()
    fixed=g.dirichlet_mask(f.X);q=(f.X-g.lower)/g.lengths
    perturb=.02*cfg['gamma']*field.length*np.prod(np.sin(np.pi*q),axis=1)
    initial=target+perturb[:,None]*np.array([1.,-.4,.7]);initial[fixed]=target[fixed]
    np.savez_compressed(OUT/f'SHEAR_{kind}_SOLVE_INPUT.npz',initial=initial,external=external,fixed=fixed)
    print('START SOLVE',kind,3*int((~fixed).sum()),flush=True)
    x,info=solve(f,material(),target,fixed,initial,external,h,sc['F_ref_N'],cfg)
    st=evaluate(f,x,material());res=st.force+external;scale=max(float(np.linalg.norm(external[~fixed])),sc['F_ref_N'])
    fr=float(np.linalg.norm(res[~fixed])/scale);e=np.linalg.norm(st.Fbar-Ftarget,axis=(1,2))/cfg['gradient_scale']
    ea=np.linalg.norm(st.Fbar-Fan,axis=(1,2))/cfg['gradient_scale'];ue=np.linalg.norm(x-target,axis=1)
    targetE=load('SHEAR_NODEWISE_IDENTITIES.json')['target_energy_J']
    result=dict(kind=kind,solve=info,force_residual=fr,solver_force_scale_N=scale,F_recovery_error_max=float(e.max()),F_analytic_error_max=float(ea.max()),
                displacement_error_max_m=float(ue.max()),displacement_relative_L2=float(np.linalg.norm(x-target)/np.linalg.norm(target-f.X)),
                initial_perturbation_max_m=normmax(initial-target),initial_x_perturbation_max_m=float(abs(initial[:,0]-target[:,0]).max()),
                energy_J=float(st.energy),target_energy_J=targetE,energy_deviation_J=float(st.energy-targetE),energy_relative_deviation=float((st.energy-targetE)/targetE),
                fixed_displacement_error_m=float(abs(x[fixed]-target[fixed]).max()),PASS=bool(info['converged'] and fr<=cfg['mechanical_residual_gate'] and (kind!='NL' or e.max()<=cfg['gradient_recovery_gate'])))
    save('SHEAR_NL_MANUFACTURED_RECOVERY.json' if kind=='NL' else 'SHEAR_CLASSICAL_MANUFACTURED_RESULT.json',result)
    np.savez_compressed(OUT/f'SHEAR_{kind}_RECOVERED_STATE.npz',X=f.X,x=x,target=target,initial=initial,volume=f.volume,external=external,
                        force_internal=st.force,reaction=-res*fixed[:,None],fixed=fixed,Fbar=st.Fbar,F_recovery_error=e,F_analytic_error=ea,u_error=ue)
    bins=np.rint(g.distance(f.X)/h).astype(int)
    csvsave(f'SHEAR_{kind}_RECOVERY_PROFILE.csv',[dict(depth_h=int(k),depth_over_delta=k*h/f.horizon,nodes=int((bins==k).sum()),
                F_recovery_max=float(e[bins==k].max()),F_recovery_RMS=float(np.sqrt(np.mean(e[bins==k]**2))),F_analytic_max=float(ea[bins==k].max()),
                u_error_max_m=float(ue[bins==k].max())) for k in np.unique(bins)])
    print('DONE',kind,'force',fr,'F',e.max(),flush=True)
    if kind=='NL':require(result['PASS'],'NL manufactured shear recovery')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True,choices=['reference_affine','oracle','decomposition','n27','NL','classical']);args=ap.parse_args()
    cfg=load('PREDECLARED_CONFIG.json');require(integrity()['PASS'],'protected input/source integrity')
    require(psutil.virtual_memory().available>=.2*psutil.virtual_memory().total,'RAM headroom before stage')
    started=time.monotonic();stop=threading.Event();process=psutil.Process();stats=dict(stage=args.stage,oracle_workers=8 if args.stage in ['reference_affine','oracle'] else 0,
        equilibrium_workers=1 if args.stage in ['NL','classical'] else 0,BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',peak_process_tree_RSS_bytes=0,minimum_available_RAM_bytes=psutil.virtual_memory().available)
    def monitor():
        while not stop.is_set():
            rss=0
            for p in [process]+process.children(recursive=True):
                try:rss+=p.memory_info().rss
                except psutil.Error:pass
            stats['peak_process_tree_RSS_bytes']=max(stats['peak_process_tree_RSS_bytes'],rss)
            stats['minimum_available_RAM_bytes']=min(stats['minimum_available_RAM_bytes'],psutil.virtual_memory().available);stop.wait(.1)
    t=threading.Thread(target=monitor,daemon=True);t.start()
    try:
        if args.stage=='reference_affine':reference_affine(cfg)
        elif args.stage=='oracle':oracle(cfg)
        elif args.stage=='decomposition':decompose(cfg)
        elif args.stage=='n27':n27(cfg)
        else:recovery(cfg,args.stage)
    finally:
        stop.set();t.join();stats['wall_seconds']=time.monotonic()-started;stats['peak_main_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        stats['minimum_available_RAM_fraction']=stats['minimum_available_RAM_bytes']/psutil.virtual_memory().total
        save(f'RESOURCES_{args.stage}.json',stats);save('HISTORICAL_INTEGRITY.json',integrity())
    require(stats['minimum_available_RAM_fraction']>=.2,'RAM headroom during stage')


if __name__=='__main__':main()

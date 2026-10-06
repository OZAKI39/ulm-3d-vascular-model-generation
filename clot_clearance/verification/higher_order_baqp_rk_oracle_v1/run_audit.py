"""Fixed-state higher-order kinematics only. No stress, Z, force or equilibrium."""
import csv,gc,hashlib,json,resource,threading,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import psutil
from pd_clot.pfpd_candidate.kernel import radial
from pd_clot.pfpd_candidate.higher_order_rk.basis import NAMES,physical_polynomials
from pd_clot.pfpd_candidate.higher_order_rk.rk_weights import directed_topology,build_weights
from pd_clot.pfpd_candidate.higher_order_rk.rk_kinematics import deformation,gradient,c1_baseline
from pd_clot.pfpd_candidate.higher_order_rk.manufactured_fields import exact
from pd_clot.pfpd_candidate.higher_order_rk.bond_diagnostic import associated_gradient,psi,metrics
from pd_clot.pfpd_candidate.higher_order_rk.oracle_compare import representative_nodes

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
OLD=ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'
def read(name):return json.loads((OUT/name).read_text())
def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v))
def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,default=plain,allow_nan=False)+'\n')
def csvsave(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(a):return hashlib.sha256(np.asarray(a).tobytes()).hexdigest()
def check(ok,reason):
    if not ok:save('HARD_STOP.json',{'reason':reason,'thresholds_unchanged':True});raise RuntimeError(reason)
def integrity():
    old=read('PROTECTED_INPUT_HASHES.json');changed=[p for p,h in old.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    return dict(protected=len(old),changed=changed,PASS=not changed)
def norm(e):return np.sqrt(np.sum(e*e,axis=tuple(range(1,e.ndim))))
def gate(m,cfg,channel='F'):
    return m['max_scaled']<=cfg[channel+'_max_gate'] and m['weighted_RMS_scaled']<=cfg[channel+'_RMS_gate']


def main_run(cfg):
    check(integrity()['PASS'],'Protected input changed before audit')
    with np.load(OLD/'MOMENT_TENSORS.npz') as a:
        X=a['X'];V=a['volume'];depth=a['distance_face_m']/cfg['delta_m'];classes=a['geometry_class'];family_size=a['family_size']
    classes=np.where(classes=='bulk-like','interior-like',classes)
    with np.load(OLD/'BOND_GEOMETRY.npz') as a:pairs=a['pairs'].astype(np.int64);bw=a['bond_weight'];bd=a['midpoint_depth_over_delta']
    L=cfg['L_geometry_m'];h=cfg['h_m'];delta=cfg['delta_m'];i,j=pairs.T;xi=X[j]-X[i];omega=radial(norm(xi)/delta)
    ident=json.loads((OLD/'FAMILY_IDENTITY.json').read_text())
    for key,a in [('X',X),('volume',V),('pairs',pairs),('omega',omega)]:check(sha(a)==ident[key+'_sha256'],'Saved family mismatch: '+key)
    source,target,ptr=directed_topology(pairs,len(X));dw=radial(norm(X[target]-X[source])/delta)
    labels=['interior-like','face-near','edge-near','corner-near'];node_bins=np.rint(depth*delta/h).astype(int);bond_bins=np.rint(bd*delta/h*2).astype(int)
    np.savez_compressed(OUT/'REFERENCE_GEOMETRY.npz',X=X,volume=V,pairs=pairs,omega=omega,bond_weight=bw,depth_over_delta=depth,
        bond_depth_over_delta=bd,geometry_class=classes,source=source.astype(np.int32),target=target.astype(np.int32),ptr=ptr)
    save('INPUT_IDENTITY.json',dict(previous_family=ident,all_hashes_match=True,configuration_sha256=hashlib.sha256((OUT/'PD2_CONFIGURATION.json').read_bytes()).hexdigest()))
    models={};conditioning={}
    for model,dim in [('PD2',9),('RK2_DIAG',6)]:
        data=build_weights(X,V,source,target,ptr,dw,delta,dim,cfg['CPU_workers'],cfg['block_nodes'])
        np.savez_compressed(OUT/(model+'_REFERENCE_CACHE.npz'),**data)
        cond=data['condition'];rank=data['rank'];sv=data['singular_values']
        summ=dict(min_rank=int(rank.min()),rank_deficient_count=int((rank<dim).sum()),sigma_min_min=float(sv[:,-1].min()),
            sigma_max_max=float(sv[:,0].max()),condition_median=float(np.median(cond)),condition_P95=float(np.percentile(cond,95)),
            condition_P99=float(np.percentile(cond,99)),condition_max=float(cond.max()),moment_symmetry_relative=float(norm(data['M']-data['M'].transpose(0,2,1)).max()/norm(data['M']).max()))
        conditioning[model]=summ
        if model=='PD2':
            rows=[]
            for k in range(len(X)):
                rows.append(dict(node=k,depth_over_delta=depth[k],geometry_class=classes[k],family_size=int(family_size[k]),rank=int(rank[k]),
                    sigma_min=sv[k,-1],sigma_max=sv[k,0],condition=cond[k],flag='WELL_CONDITIONED' if cond[k]<=1e6 else ('ELEVATED' if cond[k]<=1e10 else 'SEVERE'),
                    **{'sigma_'+str(c+1):sv[k,c] for c in range(dim)}))
            csvsave('PD2_NODE_CONDITIONING.csv',rows)
            grouped=[]
            for label,mask in [('all',np.ones(len(X),bool))]+[(q,classes==q) for q in labels]+[('depth_h_'+str(k),node_bins==k) for k in np.unique(node_bins)]:
                cc=cond[mask];grouped.append(dict(group=label,nodes=int(mask.sum()),rank_min=int(rank[mask].min()),rank_deficient=int((rank[mask]<dim).sum()),
                    sigma_min=float(sv[mask,-1].min()),condition_median=float(np.median(cc)),condition_P95=float(np.percentile(cc,95)),condition_P99=float(np.percentile(cc,99)),condition_max=float(cc.max())))
            csvsave('PD2_CONDITIONING_GROUPS.csv',grouped)
        save('CONDITIONING_SUMMARY.json',conditioning)
        check(summ['rank_deficient_count']==0,model+' rank deficient; preserve diagnostics and investigate unisolvency before fallback')
        models[model]=data['g'];del data;gc.collect()
    print('MOMENTS',conditioning,flush=True)
    # Affine gates precede any quadratic reconstruction.
    affine={};aff_arrays={}
    for name in ['affine_axial','affine_shear']:
        x,F=exact(name,X,L);got=deformation(X,x,source,target,models['PD2']);met=metrics(got-F,V,.004)
        affine[name]=met;aff_arrays[name]=got;aff_arrays[name+'_exact']=F
        check(met['max_scaled']<=cfg['F_max_gate'],'PD2 affine reproduction failed: '+name)
    save('PD2_AFFINE_REGRESSION.json',affine);np.savez_compressed(OUT/'PD2_AFFINE_REGRESSION.npz',**aff_arrays)
    p,pg=physical_polynomials(X);result=gradient(p,source,target,models['PD2']);residual=result-pg;scales=np.r_[np.ones(3),np.full(6,L)]
    polynorm=np.linalg.norm(residual,axis=2)/scales[None,:];poly={};prows=[]
    for a,name in enumerate(NAMES):
        mm=metrics(residual[:,a],V,scales[a]);poly[name]=mm
        check(mm['max_scaled']<=cfg['polynomial_max_gate'] and mm['weighted_RMS_scaled']<=cfg['polynomial_RMS_gate'],'PD2 polynomial reproduction failed: '+name)
        for k in range(len(X)):prows.append(dict(node=k,polynomial=name,geometry_class=classes[k],depth_over_delta=depth[k],normalization=scales[a],
            normalized_residual=polynorm[k,a],**{'residual_'+c:residual[k,a,z] for z,c in enumerate('xyz')}))
    csvsave('PD2_POLYNOMIAL_REPRODUCTION.csv',prows);save('POLYNOMIAL_SUMMARY.json',poly)
    np.savez_compressed(OUT/'POLYNOMIAL_REPRODUCTION.npz',actual=result,exact=pg,residual=residual,normalized_residual=polynorm,scales=scales)
    print('AFFINE AND ALL NINE POLYNOMIALS PASS',max(v['max_scaled'] for v in poly.values()),flush=True)
    # All three nodal gates must pass before constructing any PD2 bond state.
    fields={};summary={};profiles=[];class_rows=[];regression={}
    oldstates={'axial':('nonlocal_boundary_nonaffine_localization_v1','NODEWISE_RESIDUAL_DECOMPOSITION.npz','F_local'),
               'shear':('nonlocal_boundary_nonaffine_shear_localization_v1','NODEWISE_SHEAR_DECOMPOSITION.npz','F_analytic')}
    for name in ['axial','shear','cross']:
        x,F=exact(name,X,L);Fs={'C1':c1_baseline(X,x,V,pairs,omega)}
        for model,g in models.items():Fs[model]=deformation(X,x,source,target,g)
        if name in oldstates:
            dr,file,key=oldstates[name]
            with np.load(ROOT/'verification'/dr/file) as a:
                check(np.array_equal(x,a['x_target']) and np.array_equal(F,a[key]),'Saved manufactured state mismatch: '+name)
                diff=metrics(Fs['C1']-a['Fbar_target'],V,.004)
                check(diff['max_scaled']<=cfg['C1_regression_max_gate'],'C1 baseline regression failed: '+name)
                regression[name]=dict(**diff,target_bitwise_equal=True,F_exact_bitwise_equal=True,Fbar_bitwise_equal=bool(np.array_equal(Fs['C1'],a['Fbar_target'])))
        summary[name]={};nr=[]
        for model,got in Fs.items():
            err=got-F;met=metrics(err,V,.004);summary[name][model]={'Fbar':met}
            if model=='PD2':check(gate(met,cfg),'PD2 nodal F failed: '+name)
            en=norm(err)
            for k in np.unique(node_bins):
                q=node_bins==k;profiles.append(dict(field=name,model=model,sample='Fbar',depth_over_delta=k*h/delta,**metrics(err[q],V[q],.004)))
            for label in labels:
                q=classes==label;class_rows.append(dict(field=name,model=model,geometry_class=label,nodes=int(q.sum()),**metrics(err[q],V[q],.004)))
            for k in range(len(X)):nr.append(dict(node=k,model=model,geometry_class=classes[k],depth_over_delta=depth[k],F_error=en[k],F_error_scaled=en[k]/.004))
        filename={'axial':'PD2_AXIAL_QUADRATIC.csv','shear':'PD2_SHEAR_QUADRATIC.csv','cross':'PD2_CROSS_TERM.csv'}[name]
        csvsave(filename,nr);np.savez_compressed(OUT/(name.upper()+'_NODAL_STATE.npz'),X=X,x=x,F_exact=F,**Fs)
        fields[name]=(x,F,Fs)
    save('C1_BASELINE_REGRESSION.json',regression);csvsave('NODE_GEOMETRY_STATISTICS.csv',class_rows)
    print('C1 REGRESSION AND THREE NODAL PD2 F GATES PASS',flush=True)
    # Separate full-population Ftilde gates, with no Psi evaluation before they pass.
    all_bonds={};ftrows=[];psrows=[]
    slices=[slice(k,min(k+cfg['bond_block'],len(pairs))) for k in range(0,len(pairs),cfg['bond_block'])]
    for name,(x,F,Fs) in fields.items():
        def block(sl):
            pp=pairs[sl];mid=.5*(X[pp[:,0]]+X[pp[:,1]]);_,Fe=exact(name,mid,L)
            return {model:norm(associated_gradient(X,x,pp,Fs[model])-Fe) for model in ['C1','PD2']}
        with ThreadPoolExecutor(cfg['CPU_workers']) as pool:parts=list(pool.map(block,slices))
        bonds={model+'_F_error':np.concatenate([p[model] for p in parts]) for model in ['C1','PD2']}
        for model in ['C1','PD2']:
            err=bonds[model+'_F_error'];met=metrics(err,bw,.004);summary[name][model]['Ftilde']=met
            if model=='PD2':check(gate(met,cfg),'PD2 Ftilde failed: '+name)
            for k in np.unique(bond_bins):
                q=bond_bins==k;row=dict(field=name,model=model,sample='Ftilde',depth_over_delta=k/2*h/delta,bonds=int(q.sum()),**metrics(err[q],bw[q],.004));ftrows.append(row)
                profiles.append({key:value for key,value in row.items() if key!='bonds'})
        all_bonds[name]=bonds
    csvsave('PD2_FTILDE_DIAGNOSTIC.csv',ftrows)
    print('ALL PD2 FTILDE GATES PASS; START UNCHANGED PSI DIAGNOSTIC',flush=True)
    for name,(x,F,Fs) in fields.items():
        pe=psi(F);scale=float(pe.max());summary[name]['Psi_ref_Pa']=scale;nodepsi={}
        for model in ['C1','PD2']:
            err=psi(Fs[model])-pe;nodepsi[model+'_Psi_error']=err;met=metrics(err,V,scale);summary[name][model]['Psi_Fbar']=met
            if model=='PD2':check(gate(met,cfg,'Psi'),'PD2 nodal Psi failed: '+name)
            for k in np.unique(node_bins):
                q=node_bins==k;psrows.append(dict(field=name,model=model,sample='Fbar',depth_over_delta=k*h/delta,Psi_ref_Pa=scale,**metrics(err[q],V[q],scale)))
        np.savez_compressed(OUT/(name.upper()+'_NODE_PSI.npz'),Psi_exact=pe,**nodepsi)
        def pblock(sl):
            pp=pairs[sl];mid=.5*(X[pp[:,0]]+X[pp[:,1]]);_,Fe=exact(name,mid,L);pe=psi(Fe)
            return {model:psi(associated_gradient(X,x,pp,Fs[model]))-pe for model in ['C1','PD2']}
        with ThreadPoolExecutor(cfg['CPU_workers']) as pool:parts=list(pool.map(pblock,slices))
        for model in ['C1','PD2']:
            err=np.concatenate([p[model] for p in parts]);all_bonds[name][model+'_Psi_error']=err
            met=metrics(err,bw,scale);summary[name][model]['Psi_Ftilde']=met
            if model=='PD2':check(gate(met,cfg,'Psi'),'PD2 bond Psi failed: '+name)
            for k in np.unique(bond_bins):
                q=bond_bins==k;psrows.append(dict(field=name,model=model,sample='Ftilde',depth_over_delta=k/2*h/delta,Psi_ref_Pa=scale,**metrics(err[q],bw[q],scale)))
        np.savez_compressed(OUT/(name.upper()+'_BOND_DIAGNOSTICS.npz'),**all_bonds[name])
        if name in oldstates:
            with np.load(OLD/(name.upper()+'_BOND_DIAGNOSTICS.npz')) as b:
                d=float(np.max(abs(all_bonds[name]['C1_F_error']-b['actual_F_error']))/.004)
                ps=float(np.max(abs(all_bonds[name]['C1_Psi_error']-b['Psi_actual_bias']))/scale)
                regression[name].update(Ftilde_norm_difference_scaled=d,Psi_Ftilde_bias_difference_scaled=ps)
                check(d<=cfg['C1_regression_max_gate'] and ps<=cfg['Psi_max_gate'],'Saved C1 bond diagnostics changed: '+name)
    save('C1_BASELINE_REGRESSION.json',regression);save('C1_PD2_SUMMARY.json',summary);csvsave('F_DEPTH_PROFILES.csv',profiles);csvsave('PD2_PSI_PLUS_DIAGNOSTIC.csv',psrows)
    ids=representative_nodes(X,h);dest=OUT/'julia_crosscheck';dest.mkdir(exist_ok=True)
    for key,a in [('X',X),('volume',V),('representative_ids',ids+1),('parameters',np.array([delta,1000.,666.6666666666666]))]:np.savetxt(dest/('input_'+key+'.csv'),a,delimiter=',',fmt='%.17g')
    for name,(x,F,Fs) in fields.items():
        np.savetxt(dest/('input_'+name+'_x.csv'),x,delimiter=',',fmt='%.17g')
        np.savetxt(dest/('python_'+name+'_F.csv'),Fs['PD2'][ids].reshape(len(ids),9),delimiter=',',fmt='%.17g')
    selected=np.isin(source,ids)
    np.savez_compressed(dest/'PYTHON_REPRESENTATIVES.npz',ids=ids,source=source[selected],target=target[selected],g=models['PD2'][selected],**{k:v[2]['PD2'][ids] for k,v in fields.items()})
    # All Python gates are complete; Julia comparison remains independent and pending.
    save('PYTHON_AUDIT_STATUS.json',dict(PASS=True,force='NOT_IMPLEMENTED',production='UNCHANGED',GMLS_trigger_from_Python=False))
    print('ALL PYTHON KINEMATIC AND ENERGY-DRIVER GATES PASS',flush=True)


def main():
    cfg=read('PD2_CONFIGURATION.json');start=time.monotonic();done=threading.Event();p=psutil.Process()
    stat=dict(workers=cfg['CPU_workers'],BLAS_threads=1,peak_process_tree_RSS_bytes=0,minimum_available_RAM_fraction=1.,GPU='GPU_NOT_USED_BY_DESIGN')
    def monitor():
        while not done.is_set():
            rss=0
            for q in [p]+p.children(recursive=True):
                try:rss+=q.memory_info().rss
                except psutil.Error:pass
            stat['peak_process_tree_RSS_bytes']=max(stat['peak_process_tree_RSS_bytes'],rss)
            v=psutil.virtual_memory();stat['minimum_available_RAM_fraction']=min(stat['minimum_available_RAM_fraction'],v.available/v.total);done.wait(.05)
    check(psutil.virtual_memory().available/psutil.virtual_memory().total>=.2,'Insufficient initial RAM')
    t=threading.Thread(target=monitor,daemon=True);t.start()
    try:main_run(cfg)
    finally:
        done.set();t.join();stat['wall_seconds']=time.monotonic()-start;stat['main_peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        save('RESOURCE_USAGE.json',stat);save('HISTORICAL_INTEGRITY.json',integrity())
    check(stat['minimum_available_RAM_fraction']>=.2,'RAM headroom fell below 20%')

if __name__=='__main__':main()

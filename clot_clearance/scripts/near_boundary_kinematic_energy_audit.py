"""Fixed-state causal audit; no equilibrium, boundary redevelopment or correction."""
import csv,gc,hashlib,json,multiprocessing as mp,resource,threading,time
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor
from pathlib import Path
import numpy as np
import psutil
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.kernel import Family,radial
from pd_clot.pfpd_candidate.force_stencil_audit import lattice_offsets
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import block_family
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticShear,QuadraticAxial
from pd_clot.pfpd_candidate.nonlocal_boundary.kinematic_energy_diagnostics import (
    tensor_norm,moment_chunk,quadratic_prediction,continuous_ball,lattice_reconstruction,
    psi,node_diagnostics,bond_diagnostics,weighted_metrics,continuum_inner_energy,energy_summary)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'
SHEAR=ROOT/'verification/nonlocal_boundary_nonaffine_shear_localization_v1'
AXIAL=ROOT/'verification/nonlocal_boundary_nonaffine_localization_v1'


def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v))
def save(name,v):(OUT/name).write_text(json.dumps(v,indent=2,default=plain,allow_nan=False)+'\n')
def read(name):return json.loads((OUT/name).read_text())
def csvsave(name,rows):
    with (OUT/name).open('w',newline='') as s:
        writer=csv.DictWriter(s,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def sha(a):return hashlib.sha256(a.tobytes()).hexdigest()
def check(ok,reason):
    if not ok:save('HARD_STOP.json',dict(reason=reason,physics_and_thresholds_unchanged=True));raise RuntimeError(reason)
def integrity():
    old=read('PROTECTED_INPUT_HASHES.json');bad=[p for p,h in old.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    return dict(protected_files=len(old),changed=bad,PASS=not bad)
def material():return Material(model='svk',mu=1000.,bulk=666.6666666666666)


def run(cfg):
    h=cfg['h_m'];delta=cfg['delta_m'];L=cfg['L_geometry_m'];scale=cfg['gradient_scale'];mat=material()
    check(integrity()['PASS'],'Protected files changed before audit')
    with np.load(SHEAR/'NODEWISE_SHEAR_DECOMPOSITION.npz') as a, np.load(AXIAL/'NODEWISE_RESIDUAL_DECOMPOSITION.npz') as b:
        X=a['X'].copy();V=a['volume'].copy();targets={'shear':a['x_target'].copy(),'axial':b['x_target'].copy()}
        saved_F={'shear':a['Fbar_target'].copy(),'axial':b['Fbar_target'].copy()}
        exact_F={'shear':a['F_analytic'].copy(),'axial':b['F_local'].copy()}
        check(np.array_equal(X,b['X']) and np.array_equal(V,b['volume']),'Saved physical geometry/volume mismatch')
    geometry,f=block_family(19,L,delta)
    check(np.array_equal(X,f.X) and np.array_equal(V,f.volume),'Recreated topology input differs from saved geometry')
    fields={'shear':QuadraticShear(),'axial':QuadraticAxial()};rowidx={'shear':0,'axial':2}
    for name,field in fields.items():
        check(np.array_equal(field.position(X),targets[name]),'Saved target coordinate mismatch: '+name)
        check(np.array_equal(field.response(X[:,2])['F'],exact_F[name]),'Saved analytic gradient mismatch: '+name)
    save('FAMILY_IDENTITY.json',dict(nodes=len(X),bonds=len(f.pairs),X_sha256=sha(X),volume_sha256=sha(V),pairs_sha256=sha(f.pairs),
        omega_sha256=sha(f.omega),wbar_sha256=sha(f.wbar),saved_geometry_bitwise_equal=True,
        note='Unchanged block_family recreates topology from saved X/V; both saved exact target configurations and analytic gradients are bitwise matched. No equilibrium.'))
    full=continuous_ball(delta);save('CONTINUOUS_FULL_BALL_REFERENCE.json',full)
    check(max(full['K_relative_error'],full['M3_numerical_scaled'])<=cfg['continuous_numerical_gate'],'Continuous full-ball reference inconsistent')
    offsets=lattice_offsets(h,delta);xi=offsets*h;w=radial(np.linalg.norm(xi,axis=1)/delta)
    KL=np.einsum('i,ic,ib->cb',w*h**3,xi,xi);ML=np.einsum('i,ic,id,ib->cdb',w*h**3,xi,xi,xi);invL=np.linalg.inv(KL)
    check(set(map(tuple,offsets))==set(map(tuple,-offsets)),'Full lattice inversion symmetry missing')
    lattice=dict(neighbors=len(offsets),K=KL,M3=ML,K_relative_error_vs_continuum=float(np.linalg.norm(KL-full['K_analytic'])/np.linalg.norm(full['K_analytic'])),
        M3_norm=float(np.linalg.norm(ML)),M3_scaled=float(np.linalg.norm(ML)/(np.linalg.norm(KL)*delta)),full_reconstruction={})
    full_errors={}
    with ThreadPoolExecutor(cfg['CPU_workers']) as pool:
        for name,field in fields.items():
            args=[(X[k:k+512],xi,w,h**3,invL,field,rowidx[name]) for k in range(0,len(X),512)]
            err=np.concatenate(list(pool.map(lattice_reconstruction,args)));full_errors[name]=err
            norms=tensor_norm(err);lattice['full_reconstruction'][name]=dict(max_scaled=float(norms.max()/scale),RMS_scaled=float(np.sqrt(np.average(norms**2,weights=V))/scale))
    lattice['PASS']=lattice['M3_scaled']<=cfg['full_lattice_gate'] and all(v['max_scaled']<=cfg['full_lattice_gate'] for v in lattice['full_reconstruction'].values())
    save('FULL_LATTICE_REFERENCE.json',lattice);np.savez_compressed(OUT/'FULL_LATTICE_RECONSTRUCTION.npz',X=X,**full_errors)
    check(lattice['PASS'],'Full-lattice quadratic reference inconsistent')
    print('FULL REFERENCES PASS',lattice['K_relative_error_vs_continuum'],lattice['full_reconstruction'],flush=True)
    i,j=f.pairs.T;n=len(X);m=len(i);step=(m+cfg['CPU_workers']-1)//cfg['CPU_workers'];tasks=[]
    for k in range(0,m,step):
        sl=slice(k,k+step);tasks.append((n,i[sl],j[sl],f.xi[sl],f.omega[sl],V[i[sl]],V[j[sl]],h**3))
    with ProcessPoolExecutor(cfg['CPU_workers'],mp_context=mp.get_context('spawn')) as pool:parts=list(pool.map(moment_chunk,tasks))
    K,M3,KU,MU=[sum(p[k] for p in parts) for k in range(4)];del parts,tasks;gc.collect()
    inv=np.linalg.inv(K);invU=np.linalg.inv(KU)
    G=np.einsum('ib,ibn->in',M3[:,2,2,:],inv);GU=np.einsum('ib,ibn->in',MU[:,2,2,:],invU);GL=ML[2,2,:]@invL
    ep=cfg['beta_s']/cfg['L_shear_m']*G;eu=cfg['beta_s']/cfg['L_shear_m']*GU;ef=np.broadcast_to(cfg['beta_s']/cfg['L_shear_m']*GL,ep.shape)
    support=eu-ef;volume=ep-eu;total=ep-ef
    decomposition_error=float(np.linalg.norm(support+volume-total,axis=1).max()/scale)
    Kc=np.asarray(full['K_analytic']);bulkE=cfg['beta_s']/cfg['L_shear_m']*(M3[:,2,2,:]@np.linalg.inv(Kc))
    latticeE=cfg['beta_s']/cfg['L_shear_m']*(M3[:,2,2,:]@invL)
    den=np.linalg.norm(bulkE,axis=1);denL=np.linalg.norm(latticeE,axis=1);valid=den/scale>cfg['ratio_denominator_floor_scaled']
    modifier=np.full(n,np.nan);modifier[valid]=np.linalg.norm(ep,axis=1)[valid]/den[valid]
    modifierL=np.full(n,np.nan);vL=denL/scale>cfg['ratio_denominator_floor_scaled'];modifierL[vL]=np.linalg.norm(ep,axis=1)[vL]/denL[vL]
    axis_depth=np.minimum(X-geometry.lower,geometry.upper-X);sorted_d=np.sort(axis_depth,axis=1)
    face=sorted_d[:,0];edge=np.linalg.norm(sorted_d[:,:2],axis=1);corner=np.linalg.norm(sorted_d,axis=1)
    classes=np.array(['bulk-like','face-near','edge-near','corner-near'])[np.sum(axis_depth<delta,axis=1)]
    counts=np.bincount(i,minlength=n)+np.bincount(j,minlength=n);cond=np.linalg.cond(K)
    node_bins=np.rint(face/h).astype(int)
    np.savez_compressed(OUT/'MOMENT_TENSORS.npz',X=X,volume=V,K=K,invK=inv,M3=M3,K_uniform=KU,invK_uniform=invU,M3_uniform=MU,
        K_lattice=KL,M3_lattice=ML,K_continuous=Kc,M3_continuous=full['M3_analytic'],G_physical=G,G_uniform=GU,G_lattice=GL,
        family_size=counts,omega0=f.omega0,condition_K=cond,distance_face_m=face,distance_edge_m=edge,distance_corner_m=corner,
        geometry_class=classes,K_modifier_continuous=modifier,K_modifier_lattice=modifierL,ratio_valid=valid)
    variants={'R1_full':np.asarray(ef),'R2_support_uniform':eu,'R3_physical':ep}
    variant_rows=[]
    for name,error in variants.items():
        moments=np.broadcast_to(ML,(n,3,3,3)) if name=='R1_full' else (MU if name=='R2_support_uniform' else M3)
        normM=tensor_norm(moments);normE=np.linalg.norm(error,axis=1)
        for k in range(n):variant_rows.append(dict(variant=name,node=k,depth_over_delta=face[k]/delta,geometry_class=classes[k],
            M3_norm_m6=normM[k],F_error_x=error[k,0],F_error_y=error[k,1],F_error_z=error[k,2],F_error_norm=normE[k],F_error_scaled=normE[k]/scale))
    csvsave('FAMILY_VARIANT_COMPARISON.csv',variant_rows);del variant_rows
    def stats(error):
        q=np.linalg.norm(error,axis=1);return dict(max_scaled=float(q.max()/scale),weighted_RMS_scaled=float(np.sqrt(np.average(q*q,weights=V))/scale))
    source=dict(variants={k:stats(v) for k,v in variants.items()},support_contribution=stats(support),volume_contribution=stats(volume),total_contribution=stats(total),
        decomposition_residual_scaled=decomposition_error,K_modifier_continuous_range=[float(np.nanmin(modifier)),float(np.nanmax(modifier))],
        K_modifier_lattice_range=[float(np.nanmin(modifierL)),float(np.nanmax(modifierL))],modifier_valid_count=int(valid.sum()),
        condition_K_range=[float(cond.min()),float(cond.max())],max_M3_norm_m6=float(tensor_norm(M3).max()),
        interpretation='Modifier includes shape-tensor magnitude loss as well as anisotropy; M3=0 makes the quadratic error zero for any invertible K.')
    save('SOURCE_DECOMPOSITION.json',source)
    np.savez_compressed(OUT/'SOURCE_CONTRIBUTIONS.npz',support=support,volume=volume,total=total,physical_K_error=ep,continuous_K_error=bulkE,lattice_K_error=latticeE)
    check(decomposition_error<=cfg['F_max_gate'],'Support/volume decomposition failed')
    mid=.5*(X[i]+X[j]);bn=np.abs(f.xi/f.length[:,None]);bd=geometry.distance(mid);bond_bins=np.rint(2*bd/h).astype(int)
    orientation=np.minimum(np.searchsorted(cfg['orientation_bins'],bn[:,2],side='right')-1,4)
    weights=2*f.wbar*V[i]*V[j]
    np.savez_compressed(OUT/'BOND_GEOMETRY.npz',pairs=f.pairs.astype(np.int32),abs_n=bn,r_over_delta=f.length/delta,
        midpoint_depth_over_delta=bd/delta,bond_weight=weights,orientation_bin=orientation.astype(np.int8))
    closures=[];frows=[];prows=[];psrows=[];depth_rows=[];energy_rows=[];orows=[];all_nodes={};energies={};reproduced={}
    for name,field in fields.items():
        print('FIXED STATE',name,flush=True);st=evaluate(f,targets[name],mat);row=rowidx[name]
        check(np.array_equal(st.Fbar,saved_F[name]),'Saved target Fbar not reproduced: '+name)
        kerr=float(np.linalg.norm(st.K-K)/np.linalg.norm(st.K));reproduced[name]=dict(Fbar_bitwise_equal=True,K_independent_relative_error=kerr,Ftilde_sha256=sha(st.Ftilde),target_energy_J=st.energy)
        if name=='shear':
            old=json.loads((SHEAR/'SHEAR_FIELD_BIAS_SUMMARY.json').read_text())
            check(sha(st.Ftilde)==old['samples']['Ftilde']['discrete_gradient_sha256'],'Saved shear Ftilde hash mismatch')
        H=np.zeros((3,3,3));H[row,2,2]=2*.001/field.length
        predicted=quadratic_prediction(M3,inv,H);special=np.zeros_like(predicted);special[:,row,:]=.001/field.length*G
        check(float(tensor_norm(predicted-special).max()/scale)<=cfg['F_max_gate'],'Generic/specialized quadratic identity failed')
        psi_scale=float(psi(exact_F[name],mat).max())
        nd=node_diagnostics(st.Fbar,exact_F[name],predicted,V,mat,scale,psi_scale);all_nodes[name]=nd
        fm=weighted_metrics(nd['actual_F_error'],nd['predicted_F_error'],nd['F_prediction_residual'],V,scale)
        check(fm['residual_max_scaled']<=cfg['F_max_gate'] and fm['residual_RMS_scaled']<=cfg['F_RMS_gate'],'Fbar moment identity failed: '+name)
        npm=weighted_metrics(nd['Psi_actual_bias'],nd['Psi_predicted_bias'],nd['Psi_causal_residual'],V,psi_scale)
        check(npm['residual_max_scaled']<=cfg['Psi_causal_gate'],'Nodal Psi causal reconstruction failed: '+name)
        for k in range(n):
            frows.append(dict(field=name,node=k,X_m=X[k,0],Y_m=X[k,1],Z_m=X[k,2],depth_over_delta=face[k]/delta,geometry_class=classes[k],
                family_size=int(counts[k]),omega0_m3=f.omega0[k],condition_K=cond[k],M3_norm_m6=tensor_norm(M3[k:k+1])[0],M3_zz_norm_m6=np.linalg.norm(M3[k,2,2,:]),
                actual_norm=nd['actual_F_error'][k],predicted_norm=nd['predicted_F_error'][k],residual_norm=nd['F_prediction_residual'][k],residual_scaled=nd['F_prediction_residual'][k]/scale))
        node_out=dict(nd,X=X,volume=V,depth_over_delta=face/delta,geometry_class=classes)
        np.savez_compressed(OUT/f'{name.upper()}_NODE_DIAGNOSTICS.npz',**node_out)
        del node_out
        blocks=[]
        def args():
            for k in range(0,m,cfg['bond_block_size']):
                sl=slice(k,k+cfg['bond_block_size']);ii=i[sl];jj=j[sl]
                yield (f.xi[sl],st.Ftilde[sl],field.response(mid[sl,2])['F'],nd['actual_delta_F'][ii],nd['actual_delta_F'][jj],predicted[ii],predicted[jj],
                    weights[sl],st.W0[sl],mat,scale*cfg['ratio_denominator_floor_scaled'])
        # Bounded executor map: batches of eight avoid eager materialization of every tensor block.
        gen=args()
        with ThreadPoolExecutor(cfg['CPU_workers']) as pool:
            while True:
                batch=[]
                for _ in range(cfg['CPU_workers']):
                    try:batch.append(next(gen))
                    except StopIteration:break
                if not batch:break
                blocks.extend(pool.map(bond_diagnostics,batch))
        bonds={key:np.concatenate([b[key] for b in blocks]) for key in blocks[0]};del blocks,batch;gc.collect()
        bm=weighted_metrics(bonds['actual_F_error'],bonds['predicted_F_error'],bonds['projection_residual'],weights,scale)
        chain=weighted_metrics(bonds['actual_F_error'],bonds['moment_predicted_F_error'],bonds['moment_chain_residual'],weights,scale)
        check(max(bm['residual_max_scaled'],chain['residual_max_scaled'])<=cfg['F_max_gate'] and max(bm['residual_RMS_scaled'],chain['residual_RMS_scaled'])<=cfg['F_RMS_gate'],'Ftilde projection/causal identity failed: '+name)
        bpm=weighted_metrics(bonds['Psi_actual_bias'],bonds['Psi_predicted_bias'],bonds['Psi_causal_residual'],weights,psi_scale)
        check(bpm['residual_max_scaled']<=cfg['Psi_causal_gate'],'Bond Psi causal reconstruction failed: '+name)
        np.savez_compressed(OUT/f'{name.upper()}_BOND_DIAGNOSTICS.npz',**bonds)
        for k in np.unique(bond_bins):
            mask=bond_bins==k
            pr=weighted_metrics(bonds['actual_F_error'][mask],bonds['predicted_F_error'][mask],bonds['projection_residual'][mask],weights[mask],scale)
            prows.append(dict(field=name,depth_h=k/2,depth_over_delta=k/2*h/delta,bonds=int(mask.sum()),**pr))
        for sample,data,bins,sw in [('Fbar',nd,2*node_bins,V),('Ftilde',bonds,bond_bins,weights)]:
            for k in np.unique(bins):
                mask=bins==k
                met=weighted_metrics(data['Psi_actual_bias'][mask],data['Psi_predicted_bias'][mask],data['Psi_causal_residual'][mask],sw[mask],psi_scale)
                psrows.append(dict(field=name,sample=sample,depth_h=k/2,depth_over_delta=k/2*h/delta,samples=int(mask.sum()),Psi_ref_Pa=psi_scale,**met))
        for k in np.unique(node_bins):
            mask=node_bins==k;met=weighted_metrics(nd['actual_F_error'][mask],nd['predicted_F_error'][mask],nd['F_prediction_residual'][mask],V[mask],scale)
            depth_rows.append(dict(field=name,depth_h=int(k),depth_over_delta=k*h/delta,nodes=int(mask.sum()),M3_zz_norm_max_m6=float(np.linalg.norm(M3[mask,2,2,:],axis=1).max()),**met))
        for k in range(5):
            mask=orientation==k;good=mask&np.isfinite(bonds['projection_reduction_ratio'])
            orows.append(dict(field=name,bin=k,lower=cfg['orientation_bins'][k],upper=cfg['orientation_bins'][k+1],bonds=int(mask.sum()),
                actual_F_max_scaled=float(bonds['actual_F_error'][mask].max()/scale),predicted_F_max_scaled=float(bonds['predicted_F_error'][mask].max()/scale),
                actual_F_RMS_scaled=float(np.sqrt(np.average(bonds['actual_F_error'][mask]**2,weights=weights[mask]))/scale),
                ratio_valid=int(good.sum()),reduction_ratio_weighted_mean=float(np.average(bonds['projection_reduction_ratio'][good],weights=weights[good])) if good.any() else None))
        es=energy_summary(bonds,field,mat,L);repeat=energy_summary(bonds,field,mat,L);es['deterministic']=es==repeat
        check(es['deterministic'] and max(es['Gauss32_vs64_relative'])<=cfg['energy_bookkeeping_gate'],'Energy reference/determinism failed: '+name)
        check(abs(es['intact']['U0_J']-st.energy)/st.energy<=cfg['energy_bookkeeping_gate'],'BAQP energy weight convention mismatch')
        energies[name]=es
        for k in np.unique(bond_bins):
            mask=bond_bins==k;lo=max(0.,(k/2-.25)*h);hi=min(L/2,(k/2+.25)*h)
            shell=continuum_inner_energy(field,mat,L,lo)-continuum_inner_energy(field,mat,L,hi)
            for channel,prefix,q in [('intact','U',0),('positive','Uplus',1)]:
                a=float(bonds[prefix+'0_bond_J'][mask].sum());b=float(bonds[prefix+'1_bond_J'][mask].sum())
                energy_rows.append(dict(field=name,channel=channel,depth_h=k/2,depth_over_delta=k/2*h/delta,shell_lower_m=lo,shell_upper_m=hi,
                    U0_J=a,U1_J=b,U2_J=float(shell[q]),kinematic_bias_J=a-b,weight_quadrature_bias_J=b-float(shell[q]),total_bias_J=a-float(shell[q])))
        closures.append(dict(field=name,Psi_ref_Pa=psi_scale,Fbar=fm,Ftilde=bm,Ftilde_moment_chain=chain,Psi_Fbar=npm,Psi_Ftilde=bpm,
            specialized_prediction_residual_scaled=float(tensor_norm(predicted-special).max()/scale),PASS=True))
        print('CAUSAL CHAIN PASS',name,'Fbar',fm['residual_max_scaled'],'Ftilde',bm['residual_max_scaled'],'Psi',bpm['residual_max_scaled'],flush=True)
        del st,bonds,nd;gc.collect()
    a=all_nodes['axial']['actual_delta_F'][:,2,:];s=all_nodes['shear']['actual_delta_F'][:,0,:]
    mapped=np.linalg.norm(a-s,axis=1)/scale
    predicted_difference=np.linalg.norm(all_nodes['axial']['predicted_delta_F'][:,2,:]-all_nodes['shear']['predicted_delta_F'][:,0,:],axis=1)/scale
    check(mapped.max()<=cfg['geometry_mapping_gate'],'Axial/shear mapped bias mismatch')
    csvsave('AXIAL_SHEAR_GEOMETRY_CROSSCHECK.csv',[dict(node=k,depth_over_delta=face[k]/delta,axial_row_norm_scaled=np.linalg.norm(a[k])/scale,
        shear_row_norm_scaled=np.linalg.norm(s[k])/scale,mapped_difference_scaled=mapped[k],predicted_difference_scaled=predicted_difference[k],
        **{f'axial_{c}':a[k,z] for z,c in enumerate('xyz')},**{f'shear_{c}':s[k,z] for z,c in enumerate('xyz')}) for k in range(n)])
    save('AXIAL_SHEAR_CROSSCHECK_SUMMARY.json',dict(max_mapped_difference_scaled=float(mapped.max()),weighted_RMS_difference_scaled=float(np.sqrt(np.average(mapped*mapped,weights=V))),
        predicted_row_difference_scaled=float(predicted_difference.max()),same_geometric_G=True,PASS=True,
        floating_length_note='Preserves saved axial L=0.001125 and shear L=18*h (one ulp difference); both actual target arrays matched bitwise.'))
    csvsave('FBAR_MOMENT_CLOSURE.csv',frows);csvsave('FBAR_DEPTH_PROFILES.csv',depth_rows);csvsave('FTILDE_PROJECTION_CLOSURE.csv',prows)
    csvsave('PSI_BIAS_RECONSTRUCTION.csv',psrows);csvsave('BOND_ORIENTATION_PROFILES.csv',orows);csvsave('ENERGY_DEPTH_PROFILES.csv',energy_rows)
    save('ENERGY_SOURCE_DECOMPOSITION.json',energies);save('CAUSAL_IDENTITY_SUMMARY.json',closures);save('SAVED_TARGET_REPRODUCTION.json',reproduced)
    print('ALL REQUESTED CAUSAL GATES PASS',flush=True)


def main():
    cfg=read('PREDECLARED_CONFIG.json');process=psutil.Process();start=time.monotonic();done=threading.Event()
    stats=dict(CPU_workers=cfg['CPU_workers'],BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',peak_process_tree_RSS_bytes=0,minimum_available_RAM_bytes=psutil.virtual_memory().available)
    def monitor():
        while not done.is_set():
            rss=0
            for p in [process]+process.children(recursive=True):
                try:rss+=p.memory_info().rss
                except psutil.Error:pass
            stats['peak_process_tree_RSS_bytes']=max(rss,stats['peak_process_tree_RSS_bytes']);stats['minimum_available_RAM_bytes']=min(stats['minimum_available_RAM_bytes'],psutil.virtual_memory().available)
            done.wait(.05)
    check(psutil.virtual_memory().available/psutil.virtual_memory().total>=.2,'Initial RAM headroom insufficient')
    thread=threading.Thread(target=monitor,daemon=True);thread.start()
    try:run(cfg)
    finally:
        done.set();thread.join();stats['wall_seconds']=time.monotonic()-start;stats['peak_main_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        stats['minimum_available_RAM_fraction']=stats['minimum_available_RAM_bytes']/psutil.virtual_memory().total
        save('RESOURCE_USAGE.json',stats);save('HISTORICAL_INTEGRITY.json',integrity())
    check(stats['minimum_available_RAM_fraction']>=.2,'RAM headroom below required 20%')


if __name__=='__main__':main()

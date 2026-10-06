"""Task-local affine/first-variation audit. Phase and equilibrium are disabled."""
import csv,json,os,resource,time,threading
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import psutil
from scipy.spatial import cKDTree
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.boundary_force_diagnostics import (
    cube_coordinates,make_cube,continuum,force_terms,endpoint_weight,force_stats,required_force,cumulative)
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/pfpd_baqp_finite_domain_boundary_consistency_v1'
RESULTS=ROOT/'results/pfpd_baqp_finite_domain_boundary_consistency_v1'
DELTA=.00025125;EPS=.003
CASES=[('C5',5,.00025),('C7',7,.00025),('L050',9,.0005),('L075',13,.00075)]

def save(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,allow_nan=False,default=lambda a:a.tolist() if isinstance(a,np.ndarray) else a.item())+'\n')

def table(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def material():
    cfg=json.loads((ROOT/'verification/pfpd_at1_evolution_v1/provenance/PREDECLARED_CONFIG.json').read_text())['material']
    m=Material(cfg['model'],cfg['mu'],cfg['bulk'])
    assert m.model=='svk' and m.mu==1000 and m.lame==0
    return m

def close_mask(f,mask):
    result=mask.copy();i,j=f.pairs.T
    np.logical_and.at(result,i,mask[j]);np.logical_and.at(result,j,mask[i])
    return result

def particle_rows(c,st,terms,R):
    f=c.family;d=c.distances;count=np.bincount(f.pairs.ravel(),minlength=len(f.X));cond=np.linalg.cond(st.K)
    core=c.full_support;closed2=close_mask(f,core);closed3=close_mask(f,closed2)
    rows=[]
    for a in range(len(f.X)):
        row=dict(ID=a,X_m=f.X[a,0],Y_m=f.X[a,1],Z_m=f.X[a,2],z_layer=int(c.ijk[a,2]),
                 d_xminus_m=d[a,0],d_xplus_m=d[a,1],d_yminus_m=d[a,2],d_yplus_m=d[a,3],d_zminus_m=d[a,4],d_zplus_m=d[a,5],
                 d_min_m=d[a].min(),d_loading_m=d[a,4:].min(),d_over_delta=d[a].min()/DELTA,d_over_h=d[a].min()/c.h,
                 loading_d_over_delta=d[a,4:].min()/DELTA,neighbor_count=int(count[a]),omega0=f.omega0[a],K_condition=cond[a],
                 full_support=bool(core[a]),all_neighbors_full=bool(closed2[a]),all_Z_normalization_inputs_full=bool(closed3[a]))
        for k,v in [('force',st.force),('direct',terms['direct']),('Z',terms['correction']),('required_external',-st.force)]:
            row.update({k+'_'+axis+'_N':v[a,b] for b,axis in enumerate('xyz')})
        row.update(force_norm_N=np.linalg.norm(st.force[a]),force_norm_over_R=np.linalg.norm(st.force[a])/R,
                   fz_over_R=st.force[a,2]/R)
        rows.append(row)
    return rows

def first_layer(c,st,terms,endpoints,R):
    f=c.family;i,j=f.pairs.T;answers=[];bonds=[];groups=[]
    for side in (-1,1):
        layer=1 if side==-1 else c.n-2
        pid=int(np.flatnonzero((c.ijk[:,0]==c.n//2)&(c.ijk[:,1]==c.n//2)&(c.ijk[:,2]==layer))[0])
        ids=np.flatnonzero((i==pid)|(j==pid));den=float(st.force[pid,2]);samples=[]
        for b in ids:
            sign=1 if i[b]==pid else -1;other=int(j[b] if sign==1 else i[b]);vec=sign*f.xi[b];nz=abs(vec[2]/f.length[b]);rho=f.length[b]/DELTA
            zg=int(c.ijk[other,2]);surf=zg in (0,c.n-1)
            group='loaded_surface' if surf else ('same_first_interior_layer' if zg==layer else 'deeper_specimen')
            orient='[0.95,1]' if nz>=.95 else ('[0.75,0.95)' if nz>=.75 else ('[0.5,0.75)' if nz>=.5 else '[0,0.5)'))
            length='[0,0.25)' if rho<.25 else ('[0.25,0.5)' if rho<.5 else ('[0.5,0.75)' if rho<.75 else '[0.75,1]'))
            direct=float(sign*terms['pair_direct'][b,2]);z=float(sign*terms['pair_correction'][b,2])
            row=dict(side='top' if side==1 else 'bottom',particle_ID=pid,bond_ID=int(b),neighbor_ID=other,neighbor_z_layer=zg,
                     xi_x_m=vec[0],xi_y_m=vec[1],xi_z_m=vec[2],length_m=f.length[b],length_over_delta=rho,abs_nz=nz,omega=f.omega[b],
                     omega0_self=f.omega0[pid],omega0_neighbor=f.omega0[other],wbar=f.wbar[b],direct_fz_N=direct,Z_fz_N=z,total_fz_N=direct+z,
                     neighbor_group=group,orientation_group=orient,length_group=length)
            samples.append(row)
        for key in ['neighbor_group','orientation_group','length_group']:
            for label in sorted({r[key] for r in samples}):
                rr=[r for r in samples if r[key]==label];d=sum(r['direct_fz_N'] for r in rr);z=sum(r['Z_fz_N'] for r in rr)
                groups.append(dict(side=samples[0]['side'],particle_ID=pid,grouping=key,group=label,bonds=len(rr),direct_fz_N=d,Z_fz_N=z,total_fz_N=d+z,
                                   percent_of_signed_nodal_fz=100*(d+z)/den if abs(den)>1e-15*R else None))
        d=float(terms['direct'][pid,2]);z=float(terms['correction'][pid,2]);surface=sum(r['total_fz_N'] for r in samples if r['neighbor_group']=='loaded_surface')
        answers.append(dict(side=samples[0]['side'],particle_ID=pid,nodal_fz_N=den,direct_fz_N=d,Z_fz_N=z,direct_fraction=d/den,Z_fraction=z/den,
                            direct_abs_share=abs(d)/(abs(d)+abs(z)),Z_abs_share=abs(z)/(abs(d)+abs(z)),boundary_neighbor_fraction=surface/den,
                            interior_neighbor_fraction=(den-surface)/den,endpoint_contributions={k:float(v[pid,2]) for k,v in endpoints.items()},
                            endpoint_signed_fractions={k:float(v[pid,2]/den) for k,v in endpoints.items()}))
        assert abs(sum(r['total_fz_N'] for r in samples)-den)/R<1e-12
        bonds.extend(samples)
    return answers,bonds,groups

def analyze_case(spec):
    name,n,L=spec;t=time.monotonic();c=make_cube(n,L,DELTA);f=c.family;m=material();ref=continuum(m,L,EPS);R=ref['reaction_N']
    x=f.X@ref['F'].T;x[:,2]+=EPS*L/2
    st=evaluate(f,x,m);terms=force_terms(f,st);norm=np.linalg.norm(st.force)
    recon=float(np.linalg.norm(terms['total']-st.force)/norm);assert recon<=1e-12
    i,j=f.pairs.T
    ai=force_terms(f,st,.5*f.omega/f.omega0[i])['total'];aj=force_terms(f,st,.5*f.omega/f.omega0[j])['total']
    split=float(np.linalg.norm(ai+aj-st.force)/norm);assert split<=1e-12
    loaded=(c.ijk[:,2]==0)|(c.ijk[:,2]==n-1)
    edge=force_terms(f,st,endpoint_weight(f,loaded))['total'];inside=force_terms(f,st,endpoint_weight(f,~loaded))['total']
    split_surface=float(np.linalg.norm(edge+inside-st.force)/norm);assert split_surface<=1e-12
    endpoints=dict(index_i_half=ai,index_j_half=aj,loaded_surface_endpoint=edge,nonloaded_endpoint=inside)
    core=c.full_support;bulk=force_stats(st.force,core,R);collar=force_stats(st.force,~core,R);free=~loaded
    layers=[]
    for k in range(n):
        mask=c.ijk[:,2]==k;a=st.force[mask];near=min(k,n-1-k)
        layers.append(dict(z_layer=k,z_m=float(f.X[mask,2][0]),particle_count=int(mask.sum()),layer_type={0:'loaded_surface',1:'first_interior',2:'second_interior'}.get(near,'interior') if k!=n//2 else 'middle',
                           sum_fx_N=a[:,0].sum(),sum_fy_N=a[:,1].sum(),sum_fz_N=a[:,2].sum(),L1_N=np.linalg.norm(a,axis=1).sum(),L2_N=np.linalg.norm(a),
                           max_abs_fz_N=abs(a[:,2]).max(),mean_abs_fz_N=abs(a[:,2]).mean(),layer_fz_over_R=a[:,2].sum()/R,layer_L2_over_R=np.linalg.norm(a)/R,
                           direct_fz_over_R=terms['direct'][mask,2].sum()/R,Z_fz_over_R=terms['correction'][mask,2].sum()/R))
    symmetry=max(abs(layers[k]['sum_fz_N']+layers[-k-1]['sum_fz_N']) for k in range(n))/R
    top,topq=cumulative(c,-st.force,R,1);bottom,botq=cumulative(c,-st.force,R,-1)
    first,bonds,groups=first_layer(c,st,terms,endpoints,R)
    side_rows=[]
    for a in range(2):
        for s in (-1,1):
            mask=c.ijk[:,a]==(0 if s==-1 else n-1)
            side_rows.append(dict(axis='xy'[a],side=s,count=int(mask.sum()),net_vector_N=st.force[mask].sum(axis=0),net_lateral_over_R=float(np.linalg.norm(st.force[mask,:2].sum(axis=0))/R),max_lateral_over_R=float(np.linalg.norm(st.force[mask,:2],axis=1).max()/R)))
    virtual=float(-np.sum(st.force[:,2]*(f.X[:,2]+L/2)))
    masks=dict(full_support_core=core,boundary_collar=~core,top_collar=f.X[:,2]>L/2-DELTA,bottom_collar=f.X[:,2]<-L/2+DELTA,
               side_collars=np.min(c.distances[:,:4],axis=1)<DELTA)
    support={k:{**force_stats(-st.force,mask,R),'sum_required_fz_over_R':float(-st.force[mask,2].sum()/R),'L1_required_fz_over_total':float(abs(st.force[mask,2]).sum()/abs(st.force[:,2]).sum())} for k,mask in masks.items()}
    row=dict(case=name,n=n,L_m=L,h_m=c.h,delta_m=DELTA,particles=len(f.X),bonds=len(f.pairs),epsilon=EPS,
             material=dict(model=m.model,mu_Pa=m.mu,lambda_Pa=m.lame,bulk_Pa=m.bulk),phase='DISABLED',continuum=ref,
             sum_volume_over_L3=float(f.volume.sum()/L**3),face_area_ratios=c.face_areas/L**2,
             global_affine_energy_ratio=st.energy/ref['energy_J'],global_weight_sum_ratio=float(np.sum(2*f.wbar*f.volume[i]*f.volume[j])/L**3),
             affine_Fbar_max_error=float(abs(st.Fbar-ref['F']).max()),affine_Ftilde_max_error=float(abs(st.Ftilde-ref['F']).max()),
             W_max_relative_error=float(abs(st.W0/ref['W_Pa']-1).max()),P_max_relative_error=float(abs(st.P-ref['P_Pa']).max()/abs(ref['Pzz_Pa'])),
             force_decomposition_relative_error=recon,endpoint_decomposition_relative_error=split,loaded_endpoint_decomposition_relative_error=split_surface,
             free_force_stats=force_stats(st.force,free,R),bulk=bulk,collar=collar,bulk_collar_max_ratio=bulk['max_N']/collar['max_N'] if core.any() else None,
             bulk_patch_test='PASS' if core.any() and bulk['max_over_R']<=1e-10 and bulk['max_N']<=collar['max_N']*1e-4 else ('FAIL' if core.any() else 'NOT_AVAILABLE'),
             full_support_count=int(core.sum()),full_support_fraction=float(core.mean()),two_hop_closed_count=int(close_mask(f,core).sum()),three_hop_closed_count=int(close_mask(f,close_mask(f,core)).sum()),
             top_required=topq,bottom_required=botq,layer_symmetry_error_over_R=symmetry,side_boundary=side_rows,first_layer=first,
             manufactured_force_regions=support,manufactured_balance_max_over_R=float(np.linalg.norm(st.force+required_force(st),axis=1).max()/R),
             global_force_norm_over_R=float(np.linalg.norm(st.force.sum(axis=0))/R),global_moment_norm_over_RL=float(np.linalg.norm(np.cross(x,st.force).sum(axis=0))/(R*L)),
             affine_virtual_work_J=virtual,affine_virtual_work_over_RL=virtual/ref['dU_depsilon_J'],
             wall_seconds=time.monotonic()-t,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    assert abs(row['global_affine_energy_ratio']-1)<2e-11 and row['affine_Ftilde_max_error']<2e-12
    base=RESULTS/'tables';table(base/f'{name}_particles.csv',particle_rows(c,st,terms,R));table(base/f'{name}_layers.csv',layers)
    table(base/f'{name}_cumulative.csv',top+bottom);table(base/f'{name}_first_layer_bonds.csv',bonds);table(base/f'{name}_bond_groups.csv',groups)
    folder=OUT/name;folder.mkdir(exist_ok=True)
    np.savez_compressed(folder/'AFFINE_STATE.npz',X=f.X,x=x,volume=f.volume,force=st.force,direct=terms['direct'],Z_force=terms['correction'],
                        required_external=-st.force,full_support=core,Fbar=st.Fbar,K=st.K,Z=st.Z,omega0=f.omega0)
    save(folder/'SUMMARY.json',row)
    print(name,'core',row['full_support_count'],'bulk/collar',row['bulk_collar_max_ratio'],'R/Rcont',topq['ratio_to_continuum'],flush=True)
    return row

def fd_probe(task):
    seed,kind=task;c=make_cube(13,.00075,DELTA);f=c.family;m=material();ref=continuum(m,c.L,EPS);x=f.X@ref['F'].T
    st=evaluate(f,x,m);rng=np.random.default_rng(seed)
    if kind=='core_node':
        ids=np.flatnonzero(c.full_support);pid=int(ids[np.argmax(abs(st.force[ids,2]))]);v=np.zeros_like(x);v[pid,2]=c.h
    else:
        v=np.zeros_like(x);mask=c.full_support if kind=='core_random' else ~c.full_support
        v[mask]=rng.normal(size=(mask.sum(),3));v*=c.h/np.linalg.norm(v)
    predicted=float(-np.sum(st.force*v));rows=[]
    for step in (1e-4,3e-5,1e-5):
        plus=evaluate(f,x+step*v,m).energy;minus=evaluate(f,x-step*v,m).energy;fd=(plus-minus)/(2*step)
        rows.append(dict(step=step,finite_difference_J=fd,force_work_J=predicted,relative_error=abs(fd-predicted)/max(abs(predicted),1e-8*ref['energy_J'])))
    return dict(seed=seed,kind=kind,probes=rows,passed=min(r['relative_error'] for r in rows)<2e-5,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)

def main():
    t=time.monotonic();mem=psutil.virtual_memory();estimates=[]
    for name,n,L in CASES:
        X,V,ijk,areas=cube_coordinates(n,L);pairs=cKDTree(X).query_pairs(DELTA,output_type='ndarray');B=len(pairs)
        estimates.append(dict(case=name,particles=len(X),pairs=B,pairs_bytes=pairs.nbytes,estimated_persistent_array_bytes=800*B+400*len(X),expected_peak_bytes=1600*B+128*1024**2))
    peak=max(r['expected_peak_bytes'] for r in estimates);workers=min(8,max(1,int(min(.6*mem.total,mem.available-.2*mem.total)//peak)))
    assert peak<.6*mem.total and workers*peak<.6*mem.total
    config=dict(task='PFPD_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1',cases=CASES,delta_m=DELTA,epsilon=EPS,material=vars(material()),
                phase='COMPLETELY_DISABLED',continuum='SVK lambda=0; Pzz=mu*(1+epsilon)*((1+epsilon)^2-1)',
                force_decomposition_relative_tolerance=1e-12,bulk_noise_normalized_tolerance=1e-10,bulk_collar_ratio_max=1e-4,
                dominance_abs_share=.8,traction_relative_tolerance=.01,
                depth_quantile='First layer crossing of fractions of observed signed half-reaction; report nonmonotonicity and R/Rcont separately',
                FD_first_variation='Independent centered energy differences; eight fixed direction probes; no equilibrium',
                workers=workers,affine_jobs=4,FD_jobs=8,BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',
                physical_cores=psutil.cpu_count(logical=False),logical_cores=psutil.cpu_count(),RAM_total_bytes=mem.total,
                RAM_available_start_bytes=mem.available,memory_estimates=estimates,expected_concurrent_peak_bytes=workers*peak,
                memory_cap_fraction=.6,preexisting_swap_used_bytes=psutil.swap_memory().used,
                boundary_control_gate='Only after localized boundary forcing and a passing bulk patch diagnostic; never after significant unexplained core residual.',
                phase_or_activation_or_stability_campaigns=False)
    save(OUT/'CONFIG.json',config)
    stop=threading.Event();usage={'tree_peak_RSS_bytes':0,'minimum_available_bytes':mem.available,'process_tree_peak_swap_bytes':0};startswap=psutil.swap_memory()
    def monitor():
        proc=psutil.Process()
        while not stop.is_set():
            rss=swap=0
            for q in [proc]+proc.children(recursive=True):
                try:
                    rss+=q.memory_info().rss
                    status=Path(f'/proc/{q.pid}/status').read_text();line=next(l for l in status.splitlines() if l.startswith('VmSwap:'));swap+=int(line.split()[1])*1024
                except (psutil.NoSuchProcess,FileNotFoundError,StopIteration):pass
            usage['tree_peak_RSS_bytes']=max(usage['tree_peak_RSS_bytes'],rss);usage['process_tree_peak_swap_bytes']=max(usage['process_tree_peak_swap_bytes'],swap)
            usage['minimum_available_bytes']=min(usage['minimum_available_bytes'],psutil.virtual_memory().available);stop.wait(.1)
    monitor_thread=threading.Thread(target=monitor,daemon=True);monitor_thread.start()
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows=list(pool.map(analyze_case,CASES));save(OUT/'AFFINE_SUMMARY.json',rows)
            probes=list(pool.map(fd_probe,[(701+i,('core_node' if i==0 else 'core_random' if i<4 else 'collar_random')) for i in range(8)]));save(OUT/'FIRST_VARIATION.json',dict(passed=all(p['passed'] for p in probes),probes=probes))
    finally:
        stop.set();monitor_thread.join();end=psutil.swap_memory()
        usage.update(wall_seconds=time.monotonic()-t,workers=workers,affine_active_jobs=min(workers,4),first_variation_active_jobs=min(workers,8),
                     swap_in_delta_bytes=end.sin-startswap.sin,swap_out_delta_bytes=end.sout-startswap.sout,GPU='GPU_NOT_USED_BY_DESIGN')
        save(OUT/'RESOURCES.json',usage)
    assert all(p['passed'] for p in probes)

if __name__=='__main__':main()

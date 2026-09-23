"""Read every existing state; sample immutable FEM fields. NEVER integrate motion."""
import argparse,json,sys,os,time,hashlib,socket,platform,subprocess,csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from audit_math import local_tensors,candidate,dimensionless,validity,force_ratio,wall_coefficient,align_saved,stats

ENV={}
BRANCH=np.array([92.1310784,48.2637661,112.1937748])*1e-6
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(np.integer,np.bool_)):return x.item()
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    return x
def dump(p,x):Path(p).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def init(cfg):
    sys.path.insert(0,cfg['source'])
    from particle_3d.audit import read_frozen
    from particle_3d.field import FrozenFEMField
    from particle_3d.wall_geometry import WallGeometry
    from particle_3d.hydrodynamic_resistance import viscosity_from_frozen
    prov,mesh,flow,b=read_frozen(cfg['old_fem'])
    mu,fluid=viscosity_from_frozen(cfg['old_fem'])
    ENV.update(cfg=cfg,field=FrozenFEMField.from_grids(mesh,flow),mesh=mesh,boundaries=b,
               wall=WallGeometry.from_frozen(cfg['old_fem']),mu=mu,rho=fluid['density_kg_m3'],nu=fluid['kinematic_viscosity_m2_s'],fluid=fluid)
    ENV['ports']={k:np.asarray(v.points,float).mean(axis=0) for k,v in b.items() if k!='WALL'}
    bu=ENV['field'].sample(BRANCH).velocity_m_s
    ENV['prebranch']=BRANCH-bu/np.linalg.norm(bu)*12e-6 if np.isfinite(bu).all() else BRANCH+[0,0,12e-6]

def wall_sample(points):
    """Exact nearest finite-triangle geometry for a sphere; analytic narrow phase."""
    import vtk
    from particle_3d.convex_triangle import triangle_closest_many
    wall=ENV['wall'];normals=np.empty_like(points);distance=np.empty(len(points));ids=np.empty(len(points),int)
    closest=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);dist2=vtk.mutable(0.)
    for k,p in enumerate(points):
        wall._locator.FindClosestPoint(p,closest,cell,sub,dist2)
        i=int(cell);q,_=triangle_closest_many(p,wall.triangles[i:i+1]);delta=p-q[0];d=np.linalg.norm(delta)
        distance[k]=d;normals[k]=delta/d if d else wall.normal_in[i];ids[k]=i
    return distance,normals,ids

def region_masks(x,h,a,state):
    return dict(far_field=state==0,near_wall=(state>=1)&(state<=3),handoff=state==3,
                branch=np.linalg.norm(x-BRANCH,axis=1)<=10e-6,
                prebranch=np.linalg.norm(x-ENV['prebranch'],axis=1)<=6e-6,
                trunk_center=(x[:,2]>125e-6)&(h/a>=1),
                **{k.lower()+'_vicinity':np.linalg.norm(x-p,axis=1)<=8e-6 for k,p in ENV['ports'].items()})

def process(record):
    cfg=ENV['cfg'];p=Path(cfg['input_base'])/record['remote_path'] if 'input_base' in cfg else Path(record['local_path'])
    before=sha(p);assert before==record['sha256']
    meta=json.loads(p.with_suffix('.json').read_text());assert meta['velocity_role'].startswith('HELD_ON_PRECEDING_ACCEPTED_INTERVAL')
    with np.load(p,allow_pickle=False) as f:s=f['samples'].copy()
    if cfg.get('max_states'):s=s[:cfg['max_states']]
    out=Path(cfg['output'])/'states'/record['dataset'];out.mkdir(parents=True,exist_ok=True)
    if not len(s):return dict(dataset=record['dataset'],id=record['id'],count=0)
    x=align_saved(s);a=record['radius_m'];n=len(s);v=s[:,4:7]
    # Only N endpoint samples needed; force-evaluation fields are shifted by one row.
    f=ENV['field'].sample_many(s[:,1:4]);u=f.velocity_m_s.copy();g=f.velocity_gradient_s_inv.copy();cell=f.tetra_id.copy()
    u[1:]=f.velocity_m_s[:-1];g[1:]=f.velocity_gradient_s_inv[:-1];cell[1:]=f.tetra_id[:-1]
    assert np.isfinite(u).all() and np.isfinite(g).all(),(record['dataset'],record['id'])
    e,w,curl,gamma=local_tensors(g);slip=u-v;speed=np.linalg.norm(slip,axis=1)
    raw_slip=f.velocity_m_s-v;raw_speed=np.linalg.norm(raw_slip,axis=1)
    distance,normal,triangle=wall_sample(x);h=distance-a
    saved_h=s[:,14].copy();saved_h[1:]=s[:-1,14]
    gap_error=np.max(np.abs(h-saved_h));assert gap_error<2e-17,gap_error
    state=s[:,17].astype(int).copy();state[1:]=s[:-1,17].astype(int)
    sn=np.einsum('ij,ij->i',slip,normal);st=np.linalg.norm(slip-sn[:,None]*normal,axis=1)
    mu=ENV['mu'];rho=ENV['rho'];nu=ENV['nu'];diag=6*np.pi*mu*a
    lift=candidate(a,speed,gamma,mu,rho);drag=diag*speed
    zeta=wall_coefficient(a,h,mu);vn=np.einsum('ij,ij->i',v,normal);lub=zeta*np.abs(vn)
    lub_vec=-zeta[:,None]*vn[:,None]*normal;drag_vec=diag*slip
    # Solver's 512 eps * 6 backward-error scale; numerical guard, NOT a physical force threshold.
    # Exact self-scaled unconstrained matrix condition for one sphere + one normal wall block.
    # Propagate the frozen solver's backward-error scale to a conservative forward-error budget.
    condition=1+zeta/diag
    velocity_zero=512*6*np.finfo(float).eps*condition*np.maximum(np.linalg.norm(u,axis=1),np.linalg.norm(v,axis=1))
    force_zero=diag*velocity_zero
    lift_zero=candidate(a,velocity_zero,gamma,mu,rho)
    lubrication_zero=np.maximum(force_zero,zeta*velocity_zero)
    ld,ld_status=force_ratio(lift,drag,force_zero,lift_zero)
    ll,ll_status=force_ratio(lift,lub,lubrication_zero,lift_zero)
    ll[zeta==0]=np.nan;ll_status[zeta==0]=3 # NO_ACTIVE_LUBRICATION
    dims=dimensionless(a,speed,gamma,h,nu)
    cross=np.cross(slip,curl);cn=np.linalg.norm(cross,axis=1)
    direction=np.divide(cross,cn[:,None],out=np.full(cross.shape,np.nan),where=cn[:,None]>0)
    direction[speed<=velocity_zero]=np.nan
    data=dict(time_s=s[:,0],evaluation_time_s=s[:,19],dt_s=s[:,18],radius_m=np.full(n,a),
              endpoint_x_m=s[:,1:4],evaluation_x_m=x,particle_v_m_s=v,endpoint_u_m_s=f.velocity_m_s,
              local_u_m_s=u,gradient_s_inv=g,E_s_inv=e,W_s_inv=w,curl_s_inv=curl,
              tetra_id=cell,wall_triangle_id=triangle,wall_normal=normal,h_m=h,endpoint_h_m=s[:,14],state=state,
              endpoint_state=s[:,17].astype(int),shear_s_inv=gamma,slip_m_s=slip,slip_speed_m_s=speed,
              endpoint_slip_speed_m_s=raw_speed,normal_slip_m_s=sn,tangential_slip_m_s=st,
              CANDIDATE_SAFFMAN_MAGNITUDE=lift,drag_N=drag,lubrication_N=lub,lubrication_coefficient_kg_s=zeta,
              equivalent_unbalanced_contribution_N=np.linalg.norm(drag_vec+lub_vec,axis=1),
              lift_drag_ratio=ld,lift_lubrication_ratio=ll,lift_drag_status=ld_status,lift_lubrication_status=ll_status,
              force_numerical_zero_N=force_zero,velocity_numerical_zero_m_s=velocity_zero,
              lift_numerical_zero_N=lift_zero,lubrication_numerical_zero_N=lubrication_zero,
              self_scaled_resistance_condition=condition,
              candidate_direction_only=direction,
              simple_shear_scalar_mismatch=np.abs(gamma-np.linalg.norm(curl,axis=1))/np.maximum(gamma,np.finfo(float).tiny),
              **dims)
    for k,m in region_masks(x,h,a,state).items():data['region_'+k]=m
    for vtol in [.03,.1,.3]:data['descriptive_screen_'+str(vtol)]=validity(dims['Re_p'],dims['Re_G'],dims['wall_margin'],vtol)
    data['initial_diagnostic']=np.arange(n)==0
    target=out/p.name;np.savez_compressed(target,**data)
    assert sha(p)==before
    return dict(dataset=record['dataset'],id=record['id'],count=n,file=str(target),max_gap_error_m=float(gap_error),
                completed=record['completed'],outlet=record['outlet'],end_reason=record['end_reason'],
                endpoint_outside_count=int(np.sum(~f.inside_lumen)),min_gap_m=float(h.min()),
                max_lift_N=float(np.max(lift)),sample_sha_before=before,sample_sha_after=sha(p))

METRICS=['shear_s_inv','slip_speed_m_s','endpoint_slip_speed_m_s','Re_p','Re_G','CANDIDATE_SAFFMAN_MAGNITUDE',
         'drag_N','lubrication_N','lift_drag_ratio','lift_lubrication_ratio','h_over_a','wall_margin','slip_margin','shear_margin',
         'equivalent_unbalanced_contribution_N','simple_shear_scalar_mismatch']

def aggregate(results,cfg):
    out=Path(cfg['output']);pieces=[];plot=[]
    for r in results:
        if not r['count']:continue
        with np.load(r['file']) as f:
            keys=METRICS+['dt_s','state','lift_drag_status','lift_lubrication_status','initial_diagnostic']+[k for k in f.files if k.startswith('region_') or k.startswith('descriptive_screen_')]
            d={k:f[k] for k in keys};d['dataset']=np.full(r['count'],r['dataset']);d['trajectory_id']=np.full(r['count'],r['id']);pieces.append(d)
            ids=np.unique(np.r_[np.arange(0,r['count'],max(1,r['count']//60)),np.argmax(f['CANDIDATE_SAFFMAN_MAGNITUDE']),np.argmin(f['h_m'])])
            plot.append({k:f[k][ids] for k in ['evaluation_x_m','slip_m_s','candidate_direction_only','CANDIDATE_SAFFMAN_MAGNITUDE','shear_s_inv','wall_margin','lift_drag_ratio','state']})
    all={k:np.concatenate([d[k] for d in pieces]) for k in pieces[0]}
    masks={'whole':np.ones(len(all['state']),bool),**{k[7:]:v for k,v in all.items() if k.startswith('region_')}}
    masks.update({k:all['dataset']==k for k in np.unique(all['dataset'])})
    masks.update({f'{ds}__{region}':(all['dataset']==ds)&all['region_'+region] for ds in np.unique(all['dataset']) for region in ['far_field','near_wall','handoff','branch']})
    summary={group:{key:stats(all[key][m]) for key in METRICS} for group,m in masks.items()}
    weighted={group:{key:stats(all[key][m],all['dt_s'][m]) for key in METRICS} for group,m in masks.items()}
    counts={group:int(m.sum()) for group,m in masks.items()}
    guard={group:{str(int(c)):int(np.sum(all['lift_drag_status'][m]==c)) for c in [0,1,2]} for group,m in masks.items()}
    validity_summary=dict(wall_margin_above_one_count=int(np.sum(all['wall_margin']>1)),
        wall_margin_below_one_count=int(np.sum(all['wall_margin']<1)),
        screening_fractions={str(t):float(np.mean(all['descriptive_screen_'+str(t)])) for t in [.03,.1,.3]},
        near_wall_invalid_fraction=float(np.mean(all['wall_margin'][all['region_near_wall']]<1)) if np.any(all['region_near_wall']) else None,
        certified_classic_valid_fraction=0.0 if np.all(all['wall_margin']<1) else None,
        label='NO_FINITE_THRESHOLD_PROVES_VALIDITY; WALL_MARGIN_LT_ONE_FAILS_UNBOUNDED_NECESSARY_ORDERING')
    np.savez_compressed(out/'all_scalar_states.npz',**all)
    np.savez_compressed(out/'plot_points.npz',**{k:np.concatenate([d[k] for d in plot]) for k in plot[0]})
    dump(out/'statistics.json',dict(unweighted=summary,time_weighted=weighted,counts=counts,zero_guard=guard,validity=validity_summary))
    with (out/'statistics.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['group','metric','count','median','P90','P95','P99','max'])
        for group,vals in summary.items():
            for k,v in vals.items():writer.writerow([group,k]+[v[q] for q in ['count','median','P90','P95','P99','max']])
    # Complete representatives: outlet 03, strongest candidate, median complete, slow complete; avoid duplicates.
    complete=[r for r in results if r.get('completed') and r['count']]
    selected=[]
    def add(r):
        if r and r not in selected:selected.append(r)
    for outlet in ['OUTLET_03','OUTLET_02']:
        rs=[r for r in complete if r['outlet']==outlet]
        if rs:add(max(rs,key=lambda r:r['max_lift_N']))
    if complete:
        add(sorted(complete,key=lambda r:r['count'])[len(complete)//2]);add(max(complete,key=lambda r:r['count']))
    import shutil
    rd=out/'representatives';rd.mkdir(exist_ok=True)
    for r in selected:
        name=f"{r['dataset']}_{r['id']:06d}";shutil.copyfile(r['file'],rd/(name+'.npz'))
        with np.load(r['file']) as z:
            keys=['time_s','evaluation_time_s','shear_s_inv','slip_speed_m_s','CANDIDATE_SAFFMAN_MAGNITUDE','drag_N','lubrication_N','lift_drag_ratio','lift_lubrication_ratio','h_over_a','wall_margin','state']
            np.savetxt(rd/(name+'.csv'),np.column_stack([z[k] for k in keys]),delimiter=',',header=','.join(keys),comments='')
    dump(out/'representatives.json',selected)
    return summary,validity_summary

def eulerian(cfg):
    import pyvista as pv
    from particle_3d.field import FrozenFEMField
    from particle_3d.sonovue_adapter import read_sonovue
    contract,dist,_=read_sonovue(cfg['sonovue']);diam=dist.inverse_cdf(np.array([.1,.25,.5,.75,.9]));radii=diam*.5e-6
    f0=ENV['field'];mesh=ENV['mesh'];new=pv.read(cfg['new_flow']);f1=FrozenFEMField.from_grids(mesh,new)
    x=f0.points_m[f0.tetra].mean(axis=1);vol=np.abs(np.linalg.det(np.swapaxes(f0.points_m[f0.tetra][:,1:]-f0.points_m[f0.tetra][:,:1],1,2)))/6
    ids=np.arange(len(x))
    if cfg.get('max_states'):ids=ids[:cfg['max_states']]
    out=Path(cfg['output']);result={}
    for name,field in zip(cfg['flow_ids'],[f0,f1]):
        e,w,curl,gamma=local_tensors(field.gradients_s_inv[ids]);a=radii[:,None];rg=a*a*gamma[None,:]/ENV['nu']
        pref=candidate(a,np.ones((1,len(ids))),gamma[None,:],ENV['mu'],ENV['rho'])
        np.savez_compressed(out/(name+'_eulerian.npz'),x_m=x[ids],volume_m3=vol[ids],tetra_id=ids,
              shear_s_inv=gamma,curl_s_inv=curl,radius_quantiles_m=radii,Re_G=rg,saffman_prefactor_kg_s=pref)
        grid=mesh.copy();grid.cell_data.clear();grid.point_data.clear()
        if len(ids)==grid.n_cells:
            grid.cell_data['local_shear_rate_s_inv']=gamma
            surface=grid.extract_surface();surface.save(out/(name+'_shear_surface.vtp'))
        result[name]=dict(count=len(ids),shear_statistics=stats(gamma),volume_weighted_shear_statistics=stats(gamma,vol[ids]),
            Re_G_by_radius={str(q):stats(rg[i]) for i,q in enumerate([10,25,50,75,90])},
            prefactor_by_radius={str(q):stats(pref[i]) for i,q in enumerate([10,25,50,75,90])},
            slip_dependent_lift='NOT EVALUABLE WITHOUT PARTICLE SLIP' if name=='FLOW_2P0_MMPS' else 'USE_MATCHED_TRAJECTORY_AUDIT',
            field_scaled=False,flow_sha=sha(cfg['new_flow']) if name=='FLOW_2P0_MMPS' else sha(Path(cfg['old_fem'])/'frozen_reference/flow/steady_flow_stage_sv1_3q.vtu'))
    # Radius/slip/shear/wall sweeps are algebraic diagnostics, not new particle states.
    rows=[]
    for q,a in zip([10,25,50,75,90],radii):
        for slip in [0.,1e-8,1e-7,1e-6,1e-5,1e-4]:
            for gamma in [1.,10.,100.,1000.,10000.]:
                for ha in [.002,.01,.05,.1,1.,10.]:
                    h=ha*a;dim=dimensionless(a,slip,gamma,h,ENV['nu'])
                    rows.append([q,a,slip,gamma,ha,float(candidate(a,slip,gamma,ENV['mu'],ENV['rho'])),float(wall_coefficient(a,h,ENV['mu'])),float(dim['wall_margin'])])
    np.savetxt(out/'HYPOTHETICAL_ONLY_sensitivity.csv',rows,delimiter=',',header='diameter_quantile,radius_m,assumed_slip_m_s,shear_s_inv,h_over_a,CANDIDATE_SAFFMAN_MAGNITUDE,lubrication_coefficient_kg_s,wall_margin',comments='')
    dump(out/'eulerian_summary.json',dict(flows=result,radius_quantiles=dict(zip(['D10','D25','D50','D75','D90'],[dict(diameter_um=float(d),radius_m=float(a)) for d,a in zip(diam,radii)])),
         sonovue_contract=contract,sampling='ALL_CANONICAL_TETRA_CENTROIDS; P1_GRADIENT_CONSTANT_PER_CELL; BOTH_CELL_WEIGHTED_AND_VOLUME_WEIGHTED',
         branch_region=dict(center_m=BRANCH,radius_m=10e-6,prebranch_center_m=ENV['prebranch'],definition='DISPLAY_LANDMARK_ROI_DESCRIPTIVE_ONLY_NOT_ANATOMICAL_SEGMENTATION')))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('config');args=parser.parse_args();cfg=json.loads(Path(args.config).read_text())
    out=Path(cfg['output']);out.mkdir(parents=True,exist_ok=True);start=time.time()
    inventory=json.loads(Path(cfg['inventory']).read_text());records=[r for r in inventory['records'] if r['duplicate_of'] is None]
    if cfg.get('max_states'):
        records=[r for r in records if r['sample_count']>0][:1]
    commands=['hostname','nproc','lscpu','free -h','nvidia-smi']
    host=dict(hostname=socket.gethostname(),python=platform.python_version(),numpy=np.__version__,workers=cfg['workers'],gpu_used=False,
              commands={c:subprocess.run(c,shell=True,capture_output=True,text=True).stdout for c in commands})
    try:host['cpu_quota']=Path('/sys/fs/cgroup/cpu.max').read_text().strip()
    except OSError:pass
    dump(out/'compute_provenance.json',host);dump(out/'config_used.json',cfg)
    results=[]
    if cfg['workers']==1:
        init(cfg)
        for r in records:results.append(process(r))
    else:
        with ProcessPoolExecutor(max_workers=cfg['workers'],initializer=init,initargs=(cfg,)) as pool:
            for k,r in enumerate(pool.map(process,records,chunksize=4)):
                results.append(r)
                if k%100==0:print(json.dumps(dict(done=k,total=len(records),states=sum(r['count'] for r in results),seconds=time.time()-start)),flush=True)
    dump(out/'trajectory_audit_receipts.json',results);summary,valid=aggregate(results,cfg)
    if not ENV:init(cfg)
    eulerian(cfg)
    dump(out/'compute_complete.json',dict(status='PASS',elapsed_s=time.time()-start,sample_count=sum(r['count'] for r in results),
         trajectories=len(records),formal_solver_modified=False,trajectory_recomputed=False))
    print('COMPLETE',sum(r['count'] for r in results),flush=True)

if __name__=='__main__':main()

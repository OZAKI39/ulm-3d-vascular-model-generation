"""Read-only parallel campaign. Replay gate precedes all mechanism studies."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,contextlib,hashlib,json,csv,multiprocessing,resource,time,traceback,copy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from pd_clot.regularization.equilibrium import make_fixture
from pd_clot.regularization.damage_energy_audit import audit_state,decompose,thermodynamic_force,categories,counterfactuals

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'verification/damage_energy_audit_v1'
HISTORY=ROOT/'results/streaming_erosion_equilibrium_v1'
PLAN=json.loads((OUT/'provenance/PREDECLARED_PLAN.json').read_text())


def read(p):return json.loads(Path(p).read_text())


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)):return [clean(v) for v in value]
    if isinstance(value,(np.bool_,)):return bool(value)
    if isinstance(value,(np.integer,)):return int(value)
    if isinstance(value,(float,np.floating)):return float(value) if np.isfinite(value) else None
    return value


def write(path,value):
    with Path(path).open('x') as f:json.dump(clean(value),f,indent=2,allow_nan=False)


def csvfile(path,rows):
    with Path(path).open('x',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)


def cbase(small=False):
    c=read(ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json')
    if small:c['clot']['cells']=[4,3,3]
    return c


def energy(cloud,x,D,c):return audit_state(cloud,x,c['material'],c['safety'],D=D)


def save_pair(out,cl,x,D0,D1,c):
    a=energy(cl,x,D0,c);b=energy(cl,x,D1,c);p,s=decompose(cl,a,b)
    tolerance=PLAN['decomposition_absolute_tolerance_J']+PLAN['decomposition_relative_tolerance']*max(abs(a.global_energy['U_total_J']),abs(b.global_energy['U_total_J']))
    if abs(s['decomposition_residual_J'])>tolerance:raise AssertionError('Decomposition mismatch')
    np.savez_compressed(out/'particle_audit.npz',X=cl.X,x=x,fixed=cl.fixed,pairs=cl.pairs,volume=cl.volume,D_before=D0,D_after=D1,
        **{'old_'+k:v for k,v in a.particle.items()},**{'new_'+k:v for k,v in b.particle.items()},**p,**categories(cl,c,1-D0))
    rows=[dict(particle_id=i,**{k:float(v[i]) for k,v in p.items()},remaining_old=float(a.particle['remaining'][i]),remaining_new=float(b.particle['remaining'][i]),
               J_old=float(a.particle['J'][i]),J_new=float(b.particle['J'][i]),K_condition_old=float(a.particle['K_condition'][i]),K_condition_new=float(b.particle['K_condition'][i])) for i in range(len(x))]
    csvfile(out/'particles.csv',rows)
    s.update(U_before_J=a.global_energy['U_total_J'],U_after_J=b.global_energy['U_total_J'],
        U_constitutive_before_J=a.global_energy['U_constitutive_J'],U_constitutive_after_J=b.global_energy['U_constitutive_J'],
        U_stabilization_before_J=a.global_energy['U_stabilization_J'],U_stabilization_after_J=b.global_energy['U_stabilization_J'],
        minimum_support_rank=int(min(a.particle['support_rank'].min(),b.particle['support_rank'].min())),
        total_bond_damage_increment=float(np.sum(D1-D0)))
    return a,b,p,s


def fd_save(out,cl,x,D,ids,c):
    f=thermodynamic_force(cl,x,D,ids,c['material'],c['safety'],PLAN['fd_epsilons'],PLAN['fd_plateau_relative_tolerance'],PLAN['fd_plateau_absolute_tolerance_J'])
    np.savez_compressed(out/'finite_difference.npz',**f,bond_midpoints_m=(cl.X[cl.pairs[ids,0]]+cl.X[cl.pairs[ids,1]])/2)
    rows=[]
    for i,k in enumerate(ids):
        for j,e in enumerate(f['epsilons']):
            rows.append(dict(bond_id=int(k),epsilon=float(e),scheme=f['scheme'][i,j],Y_J=f['Y_J'][i,j],
                             Y_softening_J=f['Y_softening_J'][i,j],Y_reconstruction_J=f['Y_reconstruction_J'][i,j],
                             Y_stabilization_J=f['Y_stabilization_J'][i,j],F_change_per_damage=f['F_change_per_damage'][i,j],
                             stable=bool(f['stable'][i]),sign_resolved=bool(f['sign_resolved'][i]),estimate_J=f['estimate_J'][i],
                             uncertainty_J=f['uncertainty_J'][i],plateau_start_index=int(f['plateau_start'][i])))
    csvfile(out/'finite_difference.csv',rows)
    resolved=f['sign_resolved'];vals=f['estimate_J'][resolved]
    return f,dict(sampled_bonds=len(ids),stable_derivative_count=int(f['stable'].sum()),sign_resolved_count=int(resolved.sum()),
        negative_Y_count=int(f['negative_resolved'].sum()),negative_Y_fraction_of_all_sampled=float(f['negative_resolved'].mean()),
        negative_Y_fraction_of_resolved=float(f['negative_resolved'].sum()/resolved.sum()) if resolved.any() else None,
        raw_negative_fraction_at_largest_epsilon=float((f['Y_J'][:,-1]<0).mean()),
        Y_min_J=float(vals.min()) if len(vals) else None,Y_median_J=float(np.median(vals)) if len(vals) else None,Y_max_J=float(vals.max()) if len(vals) else None,
        unresolved_status='THERMODYNAMIC_FORCE_NUMERICALLY_UNRESOLVED' if not resolved.all() else 'ALL_SAMPLED_SIGNS_RESOLVED',
        diagnostic_only=True)


def affine(cl,name):
    F=np.eye(3)
    if name in ['uniaxial','nonaffine']:F=np.diag([1.01,.9975,.9975])
    elif name=='isotropic':F*=1.005
    elif name=='compression':F=np.diag([.99,1.0025,1.0025])
    elif name=='shear':F[0,1]=.01
    elif name!='undeformed':raise ValueError(name)
    x=cl.X@F.T
    if name=='nonaffine':
        i=np.argmin(np.linalg.norm(cl.X-cl.X.mean(axis=0),axis=1));x[i]+=.002*cl.spacing*np.array([1.,.3,-.2])
    return x,F


def pattern(cl,name,actual=None):
    i,j=cl.pairs.T;mid=(cl.X[i]+cl.X[j])/2;center=cl.X.mean(axis=0);pid=int(np.argmin(np.linalg.norm(cl.X-center,axis=1)))
    bid=int(np.argmin(np.linalg.norm(mid-center,axis=1)));p=np.zeros(len(i))
    if name=='uniform':p[:]=1.
    elif name=='one_bond':p[bid]=1.
    elif name=='particle_star':p[(i==pid)|(j==pid)]=1.
    elif name=='symmetric_patch':p[np.linalg.norm(mid-center,axis=1)<=1.5*cl.spacing]=1.
    elif name=='coherent_patch':
        source=np.array(cbase()['streaming']['bubble_center_m']);p[np.linalg.norm(mid-source,axis=1)<=3*cl.spacing]=1.
    elif name=='random_sparse':p[np.random.default_rng(PLAN['random_seed']).choice(len(p),max(1,len(p)//20),replace=False)]=1.
    elif name=='actual_pattern':p=actual/max(float(np.max(actual)),np.finfo(float).tiny)
    else:raise ValueError(name)
    return p


def select_bonds(cl,x,D,increment,c):
    a=energy(cl,x,D,c);b=energy(cl,x,D+increment,c);pa,_=decompose(cl,a,b)
    masks=categories(cl,c,1-D);i,j=cl.pairs.T;rng=np.random.default_rng(PLAN['random_seed'])
    def sample(ids,count=100):return np.sort(rng.choice(ids,min(count,len(ids)),replace=False))
    hotspots=np.argsort(pa['delta_U_total'])[:5]
    groups=dict(top_damage=np.argsort(-increment,kind='stable')[:100],energy_hotspot=np.flatnonzero(np.isin(i,hotspots)|np.isin(j,hotspots)),
                source_near=sample(np.flatnonzero(masks['SOURCE_NEAR'][i]|masks['SOURCE_NEAR'][j])),
                fixed_near=sample(np.flatnonzero(masks['FIXED_BASE_NEAR'][i]|masks['FIXED_BASE_NEAR'][j])),random=sample(np.arange(len(i))))
    ids=np.unique(np.concatenate(list(groups.values())))
    return ids,groups


def work(job):
    os.environ['NUMBA_CACHE_DIR']=str(OUT/f'provenance/runtime_cache/numba_{os.getpid()}')
    out=OUT/job['directory'];out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();cpu=time.process_time()
    with (out/'run.log').open('x') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        try:
            kind=job['kind'];c=cbase(kind in ['fixture','small_sensitivity']);cl=make_fixture(c)
            if kind=='replay':
                source=HISTORY/job['case']/f"damage_jump_{job['step']:03d}"
                z=np.load(source/'state.npz');old=read(source/'AUDIT.json');x=z['x'];D0=z['D_before'];D1=z['D_after']
                a,b,p,s=save_pair(out,cl,x,D0,D1,c)
                tol=PLAN['replay_absolute_tolerance_J']+PLAN['replay_relative_tolerance']*max(abs(old['U_before_damage_J']),abs(old['U_after_damage_J']))
                diff=s['delta_U_total_J']-old['delta_U_damage_J']
                write(out/'REPLAY_CHECK.json',dict(difference_J=diff,tolerance_J=tol,source=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256((source/'state.npz').read_bytes()).hexdigest()))
                if abs(diff)>tol or abs(s['U_before_J']-old['U_before_damage_J'])>tol or abs(s['U_after_J']-old['U_after_damage_J'])>tol:
                    raise AssertionError('STOP: audit implementation mismatch with saved historical energy')
                ids=np.argsort(-(D1-D0),kind='stable')[:10];fd,fs=fd_save(out,cl,x,D0,ids,c)
                s.update(case=job['case'],erosion_time=old['erosion_exposure_end_s'],macro_step=job['step'],historical_delta_U_J=old['delta_U_damage_J'],replay_difference_J=diff,
                         number_of_negative_local_Y_samples=fs['negative_Y_count'],local_Y_sample_count=len(ids),local_Y_unresolved_count=len(ids)-fs['sign_resolved_count'])
                write(out/'FD_SUMMARY.json',fs)
            elif kind=='fixture':
                x,F=affine(cl,job['deformation']);D0=np.zeros(len(cl.pairs));D1=job['amplitude']*pattern(cl,job['pattern'])
                a,b,p,s=save_pair(out,cl,x,D0,D1,c)
                s.update(deformation=job['deformation'],pattern=job['pattern'],damage_amplitude=job['amplitude'],
                         max_F_macro_error=float(np.max(abs(b.particle['F']-F))) if job['deformation']!='nonaffine' else None,
                         undeformed_high_severity=bool(job['deformation']=='undeformed' and abs(b.global_energy['U_total_J'])>PLAN['undeformed_absolute_tolerance_J']))
                if job['pattern']=='uniform':
                    scale=1-job['amplitude'];s['uniform_scaling_residual_J']=b.global_energy['U_total_J']-scale*a.global_energy['U_total_J']
            elif kind=='pattern':
                jp=HISTORY/'dt_erosion_0p01/damage_jump_001/state.npz';z=np.load(jp);x=z['x'];D0=z['D_before'];inc=z['D_after']-D0
                D1=D0+job['amplitude']*pattern(cl,job['pattern'],inc);a,b,p,s=save_pair(out,cl,x,D0,D1,c)
                uniform=D0+np.mean(D1-D0);u=energy(cl,x,uniform,c)
                s.update(pattern=job['pattern'],damage_amplitude=job['amplitude'],matched_total_damage_uniform_delta_U_J=a.global_energy['U_total_J']-u.global_energy['U_total_J'],
                         maximum_removed_integrity=float((D1-D0).max()),mean_removed_integrity=float((D1-D0).mean()))
            elif kind in ['small_sensitivity','baseline_sensitivity']:
                if kind=='small_sensitivity':
                    z=np.load(ROOT/'verification/mechanical_equilibrium_v1/adr_completion/cases/small_adr/equilibrium/state.npz');x=z['x'];D0=z['bond_D']
                    if job['state']=='affine':x,_=affine(cl,'uniaxial')
                    ids=np.arange(len(cl.pairs));groups=dict(all_bonds=ids);inc=np.zeros(len(ids))
                else:
                    step=1 if job['state']=='first' else 10;z=np.load(HISTORY/f'dt_erosion_0p01/damage_jump_{step:03d}/state.npz')
                    x=z['x'];D0=z['D_before'];inc=z['D_after']-D0;ids,groups=select_bonds(cl,x,D0,inc,c)
                f,s=fd_save(out,cl,x,D0,ids,c);write(out/'SELECTION.json',groups)
                m=categories(cl,c,1-D0)
                np.savez_compressed(out/'geometry.npz',X=cl.X,x=x,pairs=cl.pairs,fixed=cl.fixed,D=D0,actual_damage_increment=inc,**m)
                s.update(state=job['state'],total_bonds=len(cl.pairs),sampling_scope='All bonds' if kind=='small_sensitivity' else 'Stratified union; not population estimate')
            elif kind=='counterfactual':
                step=job['step'];z=np.load(HISTORY/f'dt_erosion_0p01/damage_jump_{step:03d}/state.npz');x=z['x'];D0=z['D_before'];D1=z['D_after']
                s=counterfactuals(cl,x,1-D0,1-D1,c['material'],c['safety'])
            elif kind=='alpha':
                if job['state']=='actual':
                    z=np.load(HISTORY/'dt_erosion_0p01/damage_jump_001/state.npz');x=z['x'];D0=z['D_before'];D1=z['D_after']
                else:
                    c=cbase(True);cl=make_fixture(c);x,_=affine(cl,job['state']);D0=np.zeros(len(cl.pairs));D1=.01*pattern(cl,'particle_star')
                c['material']['stabilization_alpha']*=job['factor'];a,b,p,s=save_pair(out,cl,x,D0,D1,c)
                s.update(alpha=c['material']['stabilization_alpha'],factor=job['factor'],state=job['state'],label='DIAGNOSTIC ONLY / NOT A CANDIDATE MODEL / NOT USED IN PRODUCTION')
            else:raise ValueError(kind)
            write(out/'CONFIG.json',c);write(out/'JOB.json',job)
            s.update(wall_time_s=time.perf_counter()-start,cpu_time_s=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
            write(out/'SUMMARY.json',s);print(json.dumps(clean(s),indent=2),flush=True)
            return dict(directory=job['directory'],status='COMPLETED',wall_time_s=s['wall_time_s'])
        except Exception as e:
            traceback.print_exc();write(out/'FAILED.json',dict(error=repr(e),traceback=traceback.format_exc()))
            return dict(directory=job['directory'],status='FAILED',error=repr(e))


def batch(jobs,label):
    hw=read(OUT/'provenance/HARDWARE.json');available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
    workers=max(1,min(len(jobs),hw['physical_cores'],hw['allowed_cores'],int(.8*available)//hw['worker_memory_estimate_bytes']))
    rows=[];start=time.perf_counter();print(label,'workers',workers,'jobs',len(jobs),flush=True)
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        for future in as_completed([pool.submit(work,j) for j in jobs]):
            r=future.result();rows.append(r)
            if r['status']=='FAILED':print(json.dumps(r),flush=True)
    write(OUT/f'provenance/{label}_BATCH.json',dict(workers=workers,available_memory_bytes=available,jobs=rows,wall_time_s=time.perf_counter()-start))
    if any(r['status']=='FAILED' for r in rows):raise RuntimeError('STOP: retained failed audit; inspect before any downstream stage')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['replay','mechanism']);a=parser.parse_args()
    if a.phase=='replay':
        jobs=[dict(kind='replay',case=case,step=i,directory=f'actual_run_replay/{case}/jump_{i:03d}') for case,n in [('dt_erosion_0p02',5),('dt_erosion_0p01',10),('dt_erosion_0p005',20)] for i in range(1,n+1)]
        batch(jobs,'REPLAY');rows=[read(OUT/j['directory']/'SUMMARY.json') for j in jobs]
        csvfile(OUT/'actual_run_replay/damage_jump_audit.csv',rows);write(OUT/'actual_run_replay/damage_jump_audit.json',rows)
        write(OUT/'actual_run_replay/GATE.json',dict(status='PASS',jumps=35,maximum_historical_difference_J=max(abs(r['replay_difference_J']) for r in rows)))
    else:
        assert read(OUT/'actual_run_replay/GATE.json')['status']=='PASS'
        jobs=[]
        for deformation in ['undeformed','uniaxial','isotropic','compression','shear','nonaffine']:
            for scale in PLAN['uniform_scales']:
                jobs.append(dict(kind='fixture',deformation=deformation,pattern='uniform',amplitude=1-scale,directory=f'fixtures/{deformation}/uniform_{scale:g}'))
            for p in ['one_bond','particle_star','symmetric_patch','coherent_patch','random_sparse']:
                jobs.append(dict(kind='fixture',deformation=deformation,pattern=p,amplitude=.01,directory=f'fixtures/{deformation}/{p}'))
        for p in ['uniform','one_bond','particle_star','symmetric_patch','coherent_patch','random_sparse','actual_pattern']:
            for amp in PLAN['damage_amplitudes']:jobs.append(dict(kind='pattern',pattern=p,amplitude=amp,directory=f'fixtures/pattern_study/{p}_{amp:g}'))
        for state in ['equilibrium','affine']:jobs.append(dict(kind='small_sensitivity',state=state,directory=f'bond_sensitivity/small_{state}'))
        for state in ['first','last']:jobs.append(dict(kind='baseline_sensitivity',state=state,directory=f'bond_sensitivity/baseline_{state}'))
        for step in [1,10]:jobs.append(dict(kind='counterfactual',step=step,directory=f'counterfactuals/jump_{step:03d}'))
        for state in ['actual','uniaxial','nonaffine']:
            for factor in [0,.5,1.,2.]:jobs.append(dict(kind='alpha',state=state,factor=factor,directory=f'counterfactuals/stabilization_sensitivity/{state}_{factor:g}'))
        batch(jobs,'MECHANISM')
        write(OUT/'comparisons/LOCAL_COMPUTATIONS_COMPLETE.json',dict(status='PASS',case_count=len(jobs),literature_may_start=True))


if __name__=='__main__':main()

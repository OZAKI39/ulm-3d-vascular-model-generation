"""Gated C7 activation and early stationary branch. No stability calculations."""
from pathlib import Path
from functools import lru_cache
from unittest.mock import patch
import hashlib,json,os,time,resource,sys
import numpy as np
from scripts import c5_activation_audit as old
from scripts import petsc_vi_audit as prior
from scripts import c7_reduced_vi as reduced
from pd_clot.geometry import Cloud
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.energy import Material
ROOT=old.ROOT;OUT=ROOT/'verification/pfpd_at1_c7_early_branch_v1';RESULTS=ROOT/'results/pfpd_at1_c7_early_branch_v1'
read=old.read


def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())+'\n')


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def cfg():return read(OUT/'provenance/PREDECLARED_CONFIG.json')
def arrays(path):return old.ca.read_npz(path)


@lru_cache(None)
def fixture(n):
    path=old.v2.OUT/f'geometry/C{n}.npz';a=arrays(path);X=a['X'];V=a['volume'];f=Family.from_arrays(X,V,a['pairs'],old.v2.DELTA)
    top=a['top'];bottom=X[:,2]==X[:,2].min();fp=a['face_particle'];h=old.v2.L/(n-1)
    cloud=Cloud(X,V,top|bottom,a['surface'],fp,X[fp].copy(),a['face_normal'],a['face_area'],f.pairs,f.xi,f.length,f.omega,h,old.v2.DELTA)
    m=old.ca.physics_config()['material'];return cloud,f,Material(m['model'],m['mu'],m['bulk']),top


class SavedProblem(old.Problem):
    def __init__(self,n,strain,lower=None):
        data=fixture(n)
        if lower is None:lower=np.zeros(len(data[1].pairs))
        with patch.object(old.ca,'fixture',lambda config:data):old.ca.AuditProblem.__init__(self,old.ca.physics_config(),strain,lower)
        self.n=n


def problem(n,e,lower=None):return SavedProblem(n,e,lower)


def checkpoint_key(p,x,s,previous_hash,kind):
    h=hashlib.sha256()
    for a in [x,s,p.lower]:h.update(np.asarray(a,dtype='<f8').tobytes())
    h.update(json.dumps(dict(physics=p.config,geometry=digest(old.v2.OUT/f'geometry/C{p.n}.npz'),strain=p.strain,previous=previous_hash,kind=kind,policy=cfg()),sort_keys=True).encode())
    for module in [reduced,old.spectral,old.spectral.baqp]:h.update(Path(module.__file__).read_bytes())
    return h.hexdigest()


def budget(stage,m=0):
    import psutil
    mem=psutil.virtual_memory();c=cfg();nx=735;n=nx+m;nnz=nx*n+m*(nx+1)
    arrays_est=8*(nx*nx+4*nx*m)+24*nnz
    need=c['parent_bytes']+c['workers']*c['worker_bytes']+arrays_est
    # Reclaimable current-process RSS is already included in MemAvailable's
    # subtraction, so avoid counting it twice once a pool is live.
    own=psutil.Process();rss=sum(p.memory_info().rss for p in [own]+own.children(recursive=True) if p.is_running())
    passed=mem.available+rss>=.2*mem.total+need
    row=dict(stage=stage,time=time.time(),RAM_available=mem.available,total_RAM=mem.total,reserve=.2*mem.total,estimated_arrays=arrays_est,reduced_shape=[n,n],workers=c['workers'],estimated_process_bytes=need,current_tree_RSS=rss,swap=psutil.swap_memory()._asdict(),passed=passed)
    path=OUT/'MEMORY_BUDGET.json';rows=read(path) if path.exists() else [];rows.append(row);save(path,rows)
    assert passed,'MEMORY_HEADROOM_GATE';return row


def declare():
    import psutil
    mem=psutil.virtual_memory();cores=psutil.cpu_count(logical=False);parent=768*1024**2;worker=512*1024**2
    workers=max(1,min(cores-2,int((mem.available-.2*mem.total-parent)//worker)))
    c=dict(task='PFPD_AT1_C7_EARLY_BRANCH_V1',s_tol=1e-8,q_tol=1e-7,phase_detection_tol=1e-8,q_pool_threshold=1e-3,
        pool_threshold_source=str(old.OUT/'provenance/PREDECLARED_CONFIG.json')+' : active_selection_A3_near_threshold',
        pool_neighbors='All bonds sharing at least one endpoint of seed bonds; one ring, no recursive closure',
        max_expansions=8,max_pool_fraction=.25,workers=workers,parent_bytes=parent,worker_bytes=worker,physical_cores=cores,logical_cores=psutil.cpu_count(),
        scan=[.002,.0025,.003,.0035,.004,.0045,.005],bracket_width=1e-7,max_bisections=14,
        weak_locator='After bracket, midpoint-only diagnostic bisection, max 16 further points, to obtain |min q_scaled| <= q_tol; onset interval remains the original <=1e-7 bracket',
        force_tolerance=1e-9,KKT_tolerance=1e-8,solver=prior.config()['primary'],jacobian_steps=[1e-5,3e-6,1e-6],jacobian_seed=93121,jacobian_directions=10,jacobian_tolerance=5e-5,
        intact_solver='Unchanged PETSc RSLS settings and existing mixed-strain force residual; phase frozen exactly to zero, empty phase pool. Physical gradient is NOT used as a stability Hessian.',
        first_attempts='Activation bracket upper; if failure, two successive midpoints between diagnostic weak state and failed target',
        continuation_steps=[1e-6,2e-6,5e-6],easy_iterations=4,easy_consecutive=5,max_positive_states=20,max_phase=.05,max_new_activation_fraction=.05,
        localization=dict(boundary_distance_fraction=.15,boundary_fraction=.8,localized_fraction=.01,localized_RMS_fraction=.2,partially_localized_fraction=.05),
        GPU=dict(energy_relative=1e-11,q_scaled_error=1e-10,minimum_speedup=1.3,warm_repetitions=10,sample_bonds=5000),
        stability='NOT_ASSESSED',production_changes=False)
    save(OUT/'provenance/PREDECLARED_CONFIG.json',c)
    p=problem(7,0);f=p.family;a=arrays(old.v2.OUT/'geometry/C7.npz');fa=p.cloud.face_area;norm=p.cloud.face_normal
    checks=dict(particles=len(f.X)==343,bonds=len(f.pairs)==46059,mechanical_DOF=p.nx==735,coordinates=np.array_equal(f.X,a['X']),pairs=np.array_equal(f.pairs,a['pairs']),volume=np.array_equal(f.volume,a['volume']),bounds=bool(np.allclose([f.X.min(axis=0),f.X.max(axis=0)],np.array([[-1]*3,[1]*3])*old.v2.L/2,rtol=0,atol=1e-16)),volume_sum=abs(f.volume.sum()/old.v2.L**3-1)<2e-15,top_area=abs(fa[norm[:,2]==1].sum()/old.v2.L**2-1)<2e-15,bottom_area=abs(fa[norm[:,2]==-1].sum()/old.v2.L**2-1)<2e-15)
    record=dict(particles=len(f.X),bonds=len(f.pairs),mechanical_DOF=p.nx,h_m=p.cloud.spacing,delta_m=f.horizon,delta_over_h=f.horizon/p.cloud.spacing,volume_sum_m3=float(f.volume.sum()),top_area_m2=float(fa[norm[:,2]==1].sum()),bottom_area_m2=float(fa[norm[:,2]==-1].sum()),geometry_SHA256=digest(old.v2.OUT/'geometry/C7.npz'),checks=checks,passed=all(checks.values()))
    save(OUT/'C7_GEOMETRY_CHECK.json',record);assert record['passed'],'C7_GEOMETRY_MISMATCH';budget('declare');print('DECLARED',workers,'workers',flush=True)


def classify(qs):
    return {'STRONG_LOWER':qs>cfg()['q_tol'],'WEAK_LOWER':abs(qs)<=cfg()['q_tol'],'VIOLATED_LOWER':qs < -cfg()['q_tol']}


def bound_groups(s,lower,qs,sc):
    near=abs(s-lower)<=cfg()['s_tol']
    groups={k:v&near for k,v in classify(qs).items()}
    groups['INTERIOR']=(s>lower+cfg()['s_tol'])&(s<sc-cfg()['s_tol'])
    return groups


def field_summary(p,x,s):
    Y,q,D,qs=old.q_fields(p,x,s);groups=bound_groups(s,p.lower,qs,p.sc);met=p.metrics(x,s);state=p.state(x,s)
    near=abs(s-p.lower)<=cfg()['s_tol'];interior=(s>p.lower+cfg()['s_tol'])&(s<p.sc-cfg()['s_tol']);weak=near&(abs(qs)<=cfg()['q_tol'])
    reaction=float(-state.force[p.top,2].sum())
    return dict(strain=p.strain,metrics=met,min_q=float(q.min()),q_percentiles=np.percentile(q,[0,5,50,95,100]).tolist(),min_q_scaled=float(qs.min()),q_counts={k:int(v.sum()) for k,v in groups.items()},q_IDs={k:np.flatnonzero(v).tolist() for k,v in groups.items()},reaction_N=reaction,nominal_stress_Pa=reaction/old.v2.L**2,interior_count=int(interior.sum()),weak_count=int(weak.sum()),weak_IDs=np.flatnonzero(weak).tolist(),active_count=int((s>cfg()['phase_detection_tol']).sum()),phase_bonds=p.ns),dict(x=x,s=s,lower=p.lower,Y=Y,q=q,q_scaled=qs)


def solve_pool(p,x,s,label,pool,previous_hash='',previous_weak=()):
    folder=OUT/label;key=checkpoint_key(p,x,s,previous_hash,'coupled');summary=folder/'SUMMARY.json'
    if summary.exists():
        r=read(summary);assert r['checkpoint_key']==key,'CHECKPOINT_HASH_MISMATCH';return r,arrays(folder/'STATE.npz')
    folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/'INITIAL_STATE.npz',x=x,s=s,lower=p.lower)
    start=time.perf_counter();ids=reduced.initial_pool(p,x,s,previous_weak,cfg()['q_pool_threshold']);initial_ids=ids.copy();z=p.pack(x,s);history=[];status=None
    for iteration in range(cfg()['max_expansions']+1):
        if len(ids)>.25*p.ns:status='C7_PHASE_POOL_NOT_SPARSE';break
        if p.n==7:budget(label+f'/expansion_{iteration}',len(ids))
        r,z=reduced.solve_once(p,z,ids,folder/f'pool_{iteration:02d}',pool,cfg()['workers'])
        violation=reduced.excluded_violations(p,z,ids);checked=prior.acceptance(p,z,r['reason'])
        history.append(dict(expansion=iteration,pool_size=len(ids),pool_IDs=ids.tolist(),violated_excluded_IDs=violation.tolist(),pool_reason='INITIAL_AUTOMATIC_POOL' if iteration==0 else 'EXCLUDED_KKT_VIOLATION',**r,**checked))
        save(folder/'POOL_HISTORY.json',history)
        if len(violation)==0:
            status='ACCEPTED' if checked['scientifically_accepted'] else 'STATIONARY_FAILURE';break
        expanded=reduced.expand_pool(p,ids,violation);history[-1]['added_IDs']=np.setdiff1d(expanded,ids).tolist();save(folder/'POOL_HISTORY.json',history)
        if iteration==cfg()['max_expansions']:status='MAX_POOL_EXPANSIONS';break
        ids=expanded
    x,s=p.unpack(z);row,data=field_summary(p,x,s);excluded=np.ones(p.ns,bool);excluded[ids]=False
    row.update(label=label,status=status,scientifically_accepted=status=='ACCEPTED',checkpoint_key=key,previous_hash=previous_hash,initial_pool_size=len(initial_ids),pool_size=len(ids),expansion_count=max(len(history)-1,0),pool_history=history,wall_seconds=time.perf_counter()-start,PETSc_iterations=sum(h['iterations'] for h in history),excluded_count=int(excluded.sum()),excluded_strong_count=int((excluded&(data['q_scaled']>cfg()['q_tol'])).sum()),excluded_violation_max=float(np.maximum(-data['q_scaled'][excluded],0).max()) if excluded.any() else 0.,full_KKT_bond_count=p.ns,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,worker_peak_RSS_bytes=max([h['worker_peak_RSS_bytes'] for h in history],default=0),no_oracle_initialization=True)
    data['pool']=ids;save(summary,row);np.savez_compressed(folder/'STATE.npz',**data);return row,data


def c5_regression(pool):
    rows=[]
    for label in ['first_A1','same_load_as_original']:
        target=next(t for t in prior.targets() if t['label']==label)
        x,s,e=prior.initial();p=problem(5,target['strain'],s);x=old.tracker.displacement_predictor(p,x,e)
        r,a=solve_pool(p,x,s,'c5_regression/'+label,pool,previous_hash=hashlib.sha256(x.tobytes()+s.tobytes()).hexdigest())
        reference=prior.OUT/f'continuation/PETSC_RSLS/{label}';ref=arrays(reference/'STATE.npz');rr=read(reference/'SUMMARY.json')
        errors=dict(phase_Linf=float(abs(a['s']-ref['s']).max()),scaled_state_L2=float(np.linalg.norm((a['x'][p.free]-ref['x'][p.free])/p.length_scale)),energy_relative=float(abs(r['metrics']['total_J']/rr['metrics']['total_J']-1)),reaction_relative=float(abs(r['nominal_stress_Pa']/rr['reaction_stress_Pa']-1)))
        r.update(agreement_errors=errors,reference_SHA256=digest(reference/'STATE.npz'),passed=r['scientifically_accepted'] and all(errors[k]<=v for k,v in prior.config()['agreement'].items()),automatic_positive_IDs=np.flatnonzero(a['s']>1e-8).tolist())
        save(OUT/f'c5_regression/{label}/SUMMARY.json',r);rows.append(r);save(OUT/'c5_regression/ALL.json',dict(states=rows,passed=len(rows)==2 and all(t['passed'] for t in rows)))
        print('C5 REDUCED',label,r['passed'],errors,flush=True);assert r['passed'],'C5_REDUCED_POOL_GATE_FAILED'


def intact(e,label,previous,pool):
    folder=OUT/f'intact/{label}';p=problem(7,e);s=np.zeros(p.ns);x=old.tracker.displacement_predictor(p,previous['x'],float(previous['strain']))
    previous_hash=hashlib.sha256(previous['x'].tobytes()+np.asarray(previous['strain']).tobytes()).hexdigest();key=checkpoint_key(p,x,s,previous_hash,'s0')
    if (folder/'SUMMARY.json').exists():
        row=read(folder/'SUMMARY.json');assert row['checkpoint_key']==key,'CHECKPOINT_HASH_MISMATCH';a=arrays(folder/'STATE.npz');return row,dict(a,strain=e)
    budget('intact/'+label,0);r,z=reduced.solve_once(p,p.pack(x,s),[],folder,pool,cfg()['workers']);x,s=p.unpack(z);row,a=field_summary(p,x,s)
    state=p.state(x,s);safe=True
    try:old.check_positions(old.adr.boundary_view(p.cloud,p.top,e),x,old.GUARD)
    except Exception:safe=False
    gates=dict(force=row['metrics']['force_residual']<=1e-9,fixed_boundary=bool(np.array_equal(x[p.cloud.fixed],p.boundary[p.cloud.fixed])),positive_J=bool(np.isfinite(state.Ftilde).all() and np.all(np.linalg.det(state.Ftilde)>0)),frozen_phase=bool(np.count_nonzero(s)==0),safety=safe)
    row.update(label=label,intact_accepted=all(gates.values()),gates=gates,checkpoint_key=key,previous_hash=previous_hash,**r,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    save(folder/'SUMMARY.json',row);np.savez_compressed(folder/'STATE.npz',**a,strain=e)
    print('INTACT',label,e,row['min_q'],row['metrics']['force_residual'],flush=True);assert row['intact_accepted'],'C7_INTACT_EQUILIBRIUM_FAILED';return row,dict(a,strain=e)


def activation(pool):
    source=old.OLD/'equilibrium/C7/0.00025000/STATE.npz';a=arrays(source);previous=dict(x=a['x'],strain=.00025);rows=[];low=high=None
    for k,e in enumerate(cfg()['scan']):
        r,a=intact(e,f'scan_{k:02d}',previous,pool);rows.append(r);save(OUT/'activation/SCAN.json',rows)
        if r['min_q']<0:high=(r,a);break
        low=(r,a);previous=a
    assert low and high,'C7_ACTIVATION_NOT_BRACKETED';history=[]
    for k in range(cfg()['max_bisections']):
        if high[0]['strain']-low[0]['strain']<=cfg()['bracket_width']:break
        before=[low[0]['strain'],high[0]['strain']];r,a=intact(sum(before)/2,f'bisect_{k:02d}',low[1],pool)
        if r['min_q']>=0:low=(r,a)
        else:high=(r,a)
        history.append(dict(before=before,after=[low[0]['strain'],high[0]['strain']],label=r['label']))
        save(OUT/'activation/BISECTION.json',history)
    bracket=dict(interval=[low[0]['strain'],high[0]['strain']],width=high[0]['strain']-low[0]['strain'],lower=low[0],upper=high[0],history=history,kind='INTACT_LOWER_BOUND_KKT_NOT_STABILITY')
    save(OUT/'activation/BRACKET.json',bracket);assert bracket['width']<=1e-7
    ll,hh=low,high;weak=None
    for k in range(16):
        r,a=intact((ll[0]['strain']+hh[0]['strain'])/2,f'weak_{k:02d}',ll[1],pool)
        if abs(r['min_q_scaled'])<=cfg()['q_tol']:weak=(r,a);break
        if r['min_q']>=0:ll=(r,a)
        else:hh=(r,a)
    assert weak,'C7_DIAGNOSTIC_WEAK_POINT_NOT_RESOLVED';save(OUT/'activation/WEAK.json',weak[0]);return low,high,weak


def branch(pool):
    from scripts.c7_spatial_evidence import jacobian_gate,branch_spatial
    assert read(OUT/'c5_regression/ALL.json')['passed']
    for label in ['first_A1','same_load_as_original']:assert read(OUT/f'jacobian/C5_{label}.json')['passed']
    bracket=read(OUT/'activation/BRACKET.json');weak=read(OUT/'activation/WEAK.json')
    wa=arrays(OUT/f'intact/{weak["label"]}/STATE.npz');wp=problem(7,weak['strain']);ids=reduced.initial_pool(wp,wa['x'],wa['s'],threshold=cfg()['q_pool_threshold'])
    gate=OUT/'jacobian/C7_FIRST_WEAK.json'
    if not gate.exists():jacobian_gate(wp,wa['x'],wa['s'],ids,pool,'C7_FIRST_WEAK')
    assert read(gate)['passed']
    low=bracket['lower'];lowpath=OUT/f'intact/{low["label"]}/STATE.npz';last=arrays(lowpath);last_e=low['strain'];previous_hash=digest(lowpath)
    target=bracket['upper']['strain'];rows=[];attempts=[];previous_weak=weak['q_IDs']['WEAK_LOWER'];stop=None;easy=0;step_index=0
    for attempt in range(3):
        p=problem(7,target);x=old.tracker.displacement_predictor(p,last['x'],last_e)
        r,data=solve_pool(p,x,np.zeros(p.ns),f'branch/first_{attempt:02d}',pool,previous_hash,previous_weak);attempts.append(r)
        save(OUT/'branch/FIRST_ATTEMPTS.json',attempts)
        if r['scientifically_accepted'] and r['metrics']['maximum_s']>cfg()['phase_detection_tol']:
            r['spatial']=branch_spatial(r,data);rows.append(r);last=data;last_e=target;previous_hash=digest(OUT/r['label']/'STATE.npz');previous_weak=r['weak_IDs'];break
        if r['status']=='C7_PHASE_POOL_NOT_SPARSE':stop=r['status'];break
        target=(weak['strain']+target)/2
    if not rows:stop=stop or 'C7_POSITIVE_BRANCH_NOT_ESTABLISHED'
    else:
        for k in range(1,cfg()['max_positive_states']):
            prev=rows[-1]
            if prev['metrics']['maximum_s']>=cfg()['max_phase']:stop='MAX_PHASE_0.05_REACHED';break
            if prev['PETSc_iterations']<=cfg()['easy_iterations'] and prev['expansion_count']==0:easy+=1
            else:easy=0
            if easy>=cfg()['easy_consecutive'] and step_index<len(cfg()['continuation_steps'])-1:step_index+=1;easy=0
            target=last_e+cfg()['continuation_steps'][step_index];p=problem(7,target,last['s']);x=old.tracker.displacement_predictor(p,last['x'],last_e)
            r,data=solve_pool(p,x,last['s'],f'branch/step_{k:02d}',pool,previous_hash,previous_weak)
            if not r['scientifically_accepted']:stop=r['status'];break
            r['spatial']=branch_spatial(r,data);r['load_step']=target-last_e;rows.append(r)
            if r['active_count']-prev['active_count']>cfg()['max_new_activation_fraction']*p.ns:stop='NEW_LARGE_ACTIVATION_EVENT'
            last=data;last_e=target;previous_hash=digest(OUT/r['label']/'STATE.npz');previous_weak=r['weak_IDs']
            save(OUT/'branch/ALL.json',dict(states=rows,stop_reason=stop or 'RUNNING',attempts=len(attempts)))
            print('ACCEPTED C7',len(rows),target,r['metrics']['maximum_s'],r['active_count'],r['pool_size'],flush=True)
            if stop:break
        if stop is None:stop='MAX_PHASE_0.05_REACHED' if rows[-1]['metrics']['maximum_s']>=cfg()['max_phase'] else '20_ACCEPTED_POSITIVE_STATES'
    save(OUT/'branch/ALL.json',dict(states=rows,stop_reason=stop,attempts=len(attempts),stability='NOT_ASSESSED'))


def main():
    if sys.argv[1]=='declare':declare();return
    pool=prior.vi.make_pool(cfg()['workers'])
    try:
        if sys.argv[1]=='c5':c5_regression(pool)
        elif sys.argv[1]=='activation':activation(pool)
        elif sys.argv[1]=='branch':branch(pool)
    finally:
        if pool:pool.shutdown()


if __name__=='__main__':main()

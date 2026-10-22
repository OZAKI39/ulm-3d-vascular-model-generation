"""Restricted RK45 path, exact far resistance, original near-contact fallback.

All F1/F2 are retained. RK dense-polynomial displacement bounds certify every
accepted slow interval's nonpenetration. Stage/domain failures never clip data.
"""
import time,traceback,resource
import numpy as np,h5py
from scipy.integrate import RK45
from scipy.spatial.distance import pdist
from mevo_mb.fast_sim.accelerated_water import WaterAccelerated
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from scripts.straight_tube_exp1.f10ae.schedules import schedule
from .common import *
from .field import response_context
from .metrics import endpoint
class NeedFullPhysics(Exception):pass
def safe_gap(x):return float(pdist(x).min()-2*R)
def integrate(t):
    begin=time.perf_counter();out=BASE/'lf'/t['id'];out.mkdir(parents=True,exist_ok=True)
    if (out/'SUMMARY.json').exists():
        r=read(out/'SUMMARY.json');assert r['task']==clean(t);assert sha(BASE/r['trajectory_path'])==r['trajectory_sha256'];return r
    verify();cfg=read(CONFIG)['LF'];f=t.get('frequency',2500000.);c=response_context(f);folder=BASE/'fields'/str(int(f));assert read(folder/'FIELD_GATE.json')['PASS'];basis=CachedWaterBasis(folder/('fine.h5' if t.get('variant')=='FINE' else 'coarse.h5'))
    with h5py.File(BASE/'INITIAL_STATES.h5') as h:x=h[t['state']][:]
    if t.get('offset_m'):x=x+np.asarray(t['offset_m'])
    sim=WaterAccelerated(x,c['physics'],c['response'],basis,c['cavity'],policy=policy());s=schedule(t);assert len(s.segments)==1,'LF currently admits a fixed action only';sim.control(s.segments[0].event());T=t['duration'];clock=0.;mode=False;mandatory_hf=False;records=[];times=[0.];states=[x.copy()];dense=[];nfev=0;fallback=0;transitions=[];reason=None;low_steps=0;start=x.copy();bound_min=np.inf
    maxstep=t.get('max_step_s',cfg['max_step_s']);rtol=t.get('rtol',cfg['rtol']);atol=t.get('atol_m',cfg['atol_m'])
    def clearance(z):
        # This restricted cache lies in the central bore, far from top holes.
        q=c['cavity'].lumen.local(z)
        return float(min((c['cavity'].R-np.linalg.norm(q[:,1:],axis=1)-R).min(),(c['cavity'].L-abs(q[:,0])-R).min()))
    def rhs(tt,y):
        nonlocal nfev
        z=y.reshape(x.shape);nfev+=1
        if safe_gap(z)<=.1*R:raise NeedFullPhysics('Lubrication activation reached by RK stage')
        if clearance(z)<=R:raise NeedFullPhysics('Wall safety buffer reached')
        fields=sim.fields(z)
        # Original resistance() is exactly identity beyond gap/R >= .1.
        return (fields['total_force_N']/sim.zeta).ravel()
    rk=None
    try:
        while clock<T-1e-14:
            if time.perf_counter()-begin>read(CONFIG)['budget']['per_HF_wall_s']:raise RuntimeError('LF_WALL_BUDGET')
            gap=safe_gap(x)
            if gap<=cfg['contact_guard_gap_R']*R or clearance(x)<=2*R:mode=True;rk=None
            if mode and not mandatory_hf and gap>cfg['resume_gap_R']*R and clearance(x)>4*R:mode=False
            if mode:
                if not transitions or transitions[-1]['mode']!='HF':transitions.append(dict(time_s=clock,mode='HF',minimum_gap_m=gap))
                sim.x=x.copy();h=min(cfg['HF_dt_s'],T-clock);step=sim.step(h);fallback+=1;mandatory_hf=False
                if not step['accepted']:raise RuntimeError('ORIGINAL_HF_FALLBACK_REJECTED:'+','.join(step['reasons']))
                x=sim.x.copy();clock+=h;times.append(clock);states.append(x.copy());dense.append(None);records.append(dict(time_s=clock,mode='HF',dt_s=h,min_gap_m=safe_gap(x),wall_correction_m=step['wall_correction_m'],gates=step['gates']));continue
            if not transitions or transitions[-1]['mode']!='RK45':transitions.append(dict(time_s=clock,mode='RK45',minimum_gap_m=gap))
            try:
                if rk is None:rk=RK45(rhs,clock,x.ravel(),T,rtol=rtol,atol=atol,max_step=maxstep)
                previous=x.copy();t0=clock;msg=rk.step()
                if rk.status=='failed':raise RuntimeError(msg)
                poly=rk.dense_output();Q=poly.Q.reshape(len(x),3,-1);travel=abs(poly.h)*np.linalg.norm(Q,axis=1).sum(1);bound=float(travel.max());origin_margin=np.minimum(previous-basis.origin-basis.spacing,(basis.origin+(np.array(basis.shape)-2)*basis.spacing)-previous).min()
                certified=min(safe_gap(previous)-2*bound,clearance(previous)-bound,origin_margin-bound)
                if certified<=.1*R:raise NeedFullPhysics('Conservative whole-interval geometric bound')
                x=rk.y.reshape(x.shape).copy();clock=rk.t;assert safe_gap(x)>0;wall=sim.wall_check(x)
                if wall['error_bits'] or wall['correction_m']!=0:raise RuntimeError('Original wall guard rejected RK endpoint')
                bound_min=min(bound_min,certified);times.append(clock);states.append(x.copy());dense.append(poly);low_steps+=1;records.append(dict(time_s=clock,mode='RK45',dt_s=clock-t0,min_gap_m=safe_gap(x),whole_interval_margin_m=certified,max_dense_travel_bound_m=bound))
            except NeedFullPhysics as ex:
                transitions.append(dict(time_s=clock,mode='GUARD_TO_HF',reason=str(ex)));mode=True;mandatory_hf=True;rk=None
    except Exception as ex:reason=repr(ex);(out/'exception.txt').write_text(traceback.format_exc())
    complete=reason is None and abs(clock-T)<1e-13
    # Dense output is recorded as an integrator approximation, never as raw accepted nodes.
    sample_times=np.unique(np.r_[np.linspace(0,clock,min(401,max(2,int(clock/5e-6)+1))),clock]);sampled=[]
    for tm in sample_times:
        j=min(max(np.searchsorted(times,tm,side='right')-1,0),len(times)-2)
        if len(times)==1:sampled.append(states[0]);continue
        if dense[j] is not None:sampled.append(dense[j](tm).reshape(x.shape))
        else:
            a=(tm-times[j])/(times[j+1]-times[j]);sampled.append((1-a)*states[j]+a*states[j+1])
    path=out/'trajectory.h5'
    with h5py.File(path,'w') as h:
        h['time_s']=times;h.create_dataset('positions_m',data=states,compression='gzip');h['particle_id']=np.arange(len(x));h['sample_time_s']=sample_times;h.create_dataset('sample_positions_m',data=sampled,compression='gzip');h.attrs.update(scope='LOW_FIDELITY_PREDICTION',sampling='Accepted nodes plus separately labeled SciPy RK45 dense output; HF fallback sample interpolation separately labeled',training_ready=False)
    v=sim.fields(x)['total_force_N']/sim.zeta if complete and safe_gap(x)>.1*R else None
    r=dict(id=t['id'],task=t,status='PASS' if complete else 'FAIL',reason=reason,completed_s=clock,requested_s=T,N=len(x),scope='LOW_FIDELITY_PREDICTION',trajectory_path=str(path.relative_to(BASE)),trajectory_sha256=sha(path),initial=endpoint(start,start),final=endpoint(start,x,v) if complete else None,accepted_slow_steps=low_steps,HF_fallback_steps=fallback,force_evaluations=nfev,transitions=transitions,whole_interval_min_margin_m=bound_min if np.isfinite(bound_min) else None,wall_s=time.perf_counter()-begin,rss_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,identity=identity())
    write(out/'STEP_RECORDS.json',records);write(out/'SUMMARY.json',r);print('LF',t['id'],r['status'],clock,'steps',low_steps,'HF',fallback,flush=True);return clean(r)
def batch(tasks,workers=3):
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as ex:return list(ex.map(integrate,tasks))

"""Correlated block statistics and force-density-correct periodic-flow fits."""
import math
from pathlib import Path
import numpy as np
from .common import read_json
from .units import positive


def periodic_shape(y, Ly, bin_width=0.0):
    """u = F/(m nu) * shape. F is a force on EACH particle (native x/y convention).

    Exact average of the quadratic within a bin; bins must not cross Ly/2.
    Negative force for y<Ly/2, positive for y>Ly/2.
    """
    y = np.asarray(y, float)
    h = positive(Ly,"Ly")/2
    bottom = y < h
    z = np.where(bottom, y, y-h)
    return np.where(bottom,-1,1)*(z*(h-z)-bin_width**2/12)/2


def fit_viscosity(y, velocity, force, mass, Ly, bin_width=0.0):
    positive(mass,"mass"); positive(force,"force magnitude")
    y, v = np.asarray(y,float),np.asarray(velocity,float)
    if len(y)<4 or y.shape!=v.shape or not np.isfinite(v).all():
        raise ValueError("invalid velocity samples")
    shape=periodic_shape(y,Ly,bin_width)
    X=np.column_stack((shape,np.ones(len(shape))))
    slope,offset=np.linalg.lstsq(X,v,rcond=None)[0]
    pred=X@np.array([slope,offset]); rms=float(np.sqrt(np.mean((v-pred)**2)))
    if slope<=0:
        return {"status":"INCONCLUSIVE","nu_star":None,"reason":"Nonpositive flow response", "slope":float(slope)}
    return {"status":"MEASURED","nu_star":float(force/(mass*slope)), "slope":float(slope), "offset_star":float(offset),
            "rms_star":rms,"relative_rms":rms/max(float(np.sqrt(np.mean((pred-offset)**2))),1e-30),"fit_profile_star":pred.tolist()}


def correlation_time_samples(values):
    """FFT ACF, initial positive paired sequence. tau=1 for white noise."""
    x=np.asarray(values,float)
    if len(x)<4 or not np.isfinite(x).all(): return float(max(1,len(x)))
    x=x-x.mean(); var=float(np.dot(x,x))
    if var<=1e-28: return 1.0
    size=1<<(2*len(x)-1).bit_length()
    fft=np.fft.rfft(x,n=size)
    ac=np.fft.irfft(fft*fft.conj(),n=size)[:len(x)]/var
    # Stop at the first nonpositive adjacent pair; cap by available record.
    total=0.0
    for k in range(1,len(ac)-1,2):
        pair=float(ac[k]+ac[k+1])
        if pair<=0: break
        total+=pair
    return min(float(len(x)),max(1.0,1+2*total))


def t95(df):
    # Conservative Student-t 97.5% quantiles rounded upward; no scipy dependency.
    table={1:12.707,2:4.303,3:3.183,4:2.777,5:2.571,6:2.447,7:2.365,8:2.307,9:2.263,
           10:2.229,12:2.179,15:2.132,20:2.087,30:2.043,60:2.001,120:1.980}
    return table[max(k for k in table if k<=df)] if df>=1 else None


def block_statistics(t, values, *, min_duration, min_blocks=8, fixed_block_samples=None):
    t,x=np.asarray(t,float),np.asarray(values,float)
    if len(t)<4 or x.shape!=t.shape or not np.isfinite(x).all() or np.any(np.diff(t)<=0):
        return {"status":"WINDOW_INSUFFICIENT","mean":None,"ci95_halfwidth":None,"block_means":[],"block_count":0}
    spacing=float(np.median(np.diff(t)))
    if np.max(np.abs(np.diff(t)-spacing)) > spacing*0.02:
        raise ValueError("block_statistics requires regularly sampled series")
    tau=correlation_time_samples(x)
    b=max(1,int(math.ceil(min_duration/spacing)),int(math.ceil(5*tau)),fixed_block_samples or 1)
    n=len(x)//b
    means=[float(np.mean(x[i*b:(i+1)*b])) for i in range(n)]
    centers=[float(np.mean(t[i*b:(i+1)*b])) for i in range(n)]
    sem=float(np.std(means,ddof=1)/math.sqrt(n)) if n>=2 else None
    hw=t95(n-1)*sem if sem is not None else None
    # CI is descriptive when blocks are insufficient, never promoted to PASS.
    return {"status":"SUFFICIENT" if n>=min_blocks else "WINDOW_INSUFFICIENT", "mean":float(np.mean(x)),
            "std":float(np.std(x,ddof=1)),"ci95_halfwidth":hw,"block_count":n,"block_samples":b,
            "block_duration_star":b*spacing,"block_means":means,"block_times_star":centers,
            "discarded_tail_samples":len(x)-n*b,"correlation_time_samples":tau,"effective_sample_estimate":len(x)/tau,
            "sample_count":len(x),"sample_interval_star":spacing,"interval_star":[float(t[0]),float(t[-1])],
            "ci_method":"Student t of nonoverlapping blocks; block >= 5 integrated correlation times and physical relaxation minimum"}


def native_pressure(stats, virial, mass, volume, dt=0.0):
    """Same afterIntegration data; Stats vx/vy/vz are mean MOMENTA.

    Virial plugin sums trace(stress)/3, with 0.5 per particle/pair, no V, no
    kinetic term. Compression positive. Equilibrium COM kinetic correction.
    Finite-N temperature DOF correction does not add fictitious kinetic pressure.
    """
    # Pinned stats.cu serializes currentTime after the next step increment;
    # virial_pressure.cu saves time in afterIntegration. Observed Stats time is
    # virial time + dt, although both sampled the same hook. Account for native
    # CSV %g rounding (6 significant figures), never shift the physical samples.
    def csv_quantum(t):
        t=np.abs(np.asarray(t,float))
        return np.where(t>0,10.**(np.floor(np.log10(np.where(t>0,t,1.)))-5),0.)
    if len(stats)!=len(virial):raise ValueError("PRESSURE_TIMESTAMP_MISMATCH")
    rounding=.5*(csv_quantum(stats['time'])+csv_quantum(virial['time']))+1e-10
    if not np.all(np.abs(stats['time']-(virial['time']+dt))<=rounding):
        raise ValueError("PRESSURE_TIMESTAMP_MISMATCH")
    N=stats['num_particles']; v2=sum((stats[k]/mass)**2 for k in ('vx','vy','vz'))
    k_unbiased=(stats['kBT']-mass*v2/3)*N/(N-1)
    kinetic=N/volume*(stats['kBT']-mass*v2/3)
    return virial['pressure']/volume+kinetic,k_unbiased


def read_csv(path):
    path=Path(path)
    if not path.is_file() or path.stat().st_size==0:
        raise ValueError("INPUT_MISSING "+str(path))
    a=np.genfromtxt(path,delimiter=',',names=True,dtype=float)
    a=np.atleast_1d(a)
    if a.dtype.names is None or len(a)==0 or not all(np.isfinite(a[n]).all() for n in a.dtype.names):
        raise ValueError("NONFINITE_OR_EMPTY "+str(path))
    return a


def analyze_task(directory, c, units):
    directory=Path(directory); spec=read_json(directory/'actual_parameters.json'); task=spec['task']; cand=spec['candidate']
    result={"task_id":task['id'],"candidate_id":cand['id'],"kind":task['kind'],"directory":str(directory),"parameters":spec,"status":"INCONCLUSIVE"}
    role=task.get('test_kind')
    if role is None:
        role=task['kind']
        if role=='equilibrium':
            role='eos_low' if task['n_star']<c['space']['n_star'] else ('eos_high' if task['n_star']>c['space']['n_star'] else 'equilibrium')
        elif role=='flow':
            group=[t for t in c.get('calibration',{}).get('tasks',[]) if t['candidate']==cand['id'] and t['kind']=='flow' and t['dt_factor']==1]
            if task.get('dt_factor',1)<1:role='half_dt'
            elif group and task['force_star']<max(t['force_star'] for t in group):role='half_force'
    result.update(method=cand.get('method','DPD'),test_kind=role,comparison_group=task.get('comparison_group',cand['id']))
    try:
        moments=read_csv(directory/'moments.csv'); bins=read_csv(directory/'profile_samples.csv')
        native=read_csv(directory/'native_stats.csv'); virial=read_csv(directory/'pressure/pv.csv')
        p,knative=native_pressure(native,virial,spec['m_star'],float(np.prod(spec['domain_star'])),spec['dt_star'])
    except (OSError,ValueError) as exc:
        result.update(reason=str(exc),status="INCONCLUSIVE")
        return result
    t=moments['time_star']; dt=spec['dt_star']; tol=c['proposed_acceptance']; Ly=spec['domain_star'][1]
    tau_prior=Ly**2/(4*math.pi**2*cand['nu_prior_star'])
    minimum_warm=max(tol['warmup_relaxation_multiples']*tau_prior,t[-1]*tol['minimum_warmup_fraction'])
    nb=int(round(Ly/c['space']['output_y_bin_star'])); by=bins['ux'].reshape(-1,nb)
    # Each instantaneous profile is measured, not the analytical target.
    if len(by)!=len(t):
        result['reason']='profile/moment sample count mismatch';return result
    y=bins['y_star'][:nb]
    slope_signal=(by@periodic_shape(y,Ly,c['space']['output_y_bin_star']))/np.dot(periodic_shape(y,Ly,c['space']['output_y_bin_star']),periodic_shape(y,Ly,c['space']['output_y_bin_star']))
    # Data-dependent warmup: earliest allowed cut with consistent two halves.
    chosen=None; drift=None
    signal=slope_signal if task['kind']=='flow' else moments['kBT_thermal_star']
    for frac in (0.2,0.3,0.4,0.5):
        cutoff=max(minimum_warm,t[-1]*frac)
        ids=np.flatnonzero(t>=cutoff)
        if len(ids)<16: continue
        mid=len(ids)//2
        first,last=np.mean(signal[ids[:mid]]),np.mean(signal[ids[mid:]])
        this_drift=abs(float(last-first))/max(abs(float(last)),1e-20)
        # A steady mean flow does not make a still-relaxing thermostat stationary.
        first_T,last_T=np.mean(moments['kBT_thermal_star'][ids[:mid]]),np.mean(moments['kBT_thermal_star'][ids[mid:]])
        thermal_drift=abs(float(last_T-first_T))/max(abs(float(last_T)),1e-20)
        this_drift=max(this_drift,thermal_drift)
        chosen,drift=ids,this_drift
        if this_drift<=tol['max_stationarity_drift']:break
    if chosen is None:
        result.update(status='WINDOW_INSUFFICIENT',reason='Cannot cover relaxation-based warmup and statistics', duration_star=float(t[-1]), minimum_warmup_star=minimum_warm)
        return result
    sel=chosen; start=float(t[sel[0]]); ts=t[sel]; physical_span=float((ts[-1]-ts[0])*units.t0)
    min_block=tol['block_relaxation_multiples']*tau_prior
    temp=block_statistics(ts,moments['kBT_thermal_star'][sel],min_duration=min_block,min_blocks=tol['min_blocks'])
    ps=(native['time']>=start)&(native['time']<=t[-1])
    pressure=block_statistics(native['time'][ps],p[ps],min_duration=min_block,min_blocks=tol['min_blocks'])
    result.update(status='MEASURED', actual_n_star=float(moments['N'][0]/np.prod(spec['domain_star'])),
                  actual_N=int(moments['N'][0]),mass_conserved=bool(np.all(moments['N']==moments['N'][0])),
                  density_scope='N*m/V constructed in a periodic fixed volume; no vessel density-field validation',
                  warmup_interval_star=[0,start],statistics_interval_star=[start,float(t[-1])],statistics_duration_s=physical_span,
                  relaxation_prior_star=tau_prior,stationarity_relative_drift=drift,
                  stationarity_status='PASS_PROPOSED' if drift<=tol['max_stationarity_drift'] else 'UNSTABLE',
                  temperature=temp,pressure=pressure,
                  pressure_precision_status='PASS_PROPOSED' if pressure['status']=='SUFFICIENT' and pressure['ci95_halfwidth'] is not None and pressure['ci95_halfwidth']/max(abs(pressure['mean']),1e-30)<=tol['eos_ci_relative_halfwidth'] else 'INCONCLUSIVE',
                  pressure_timestamp_note='Stats reports time at serializeAndSend, virial saves time at afterIntegration: native Stats labels are virial labels + dt. Matched same hook with explicit dt and native CSV rounding tolerance.',
                  temperature_status='PASS_PROPOSED' if temp['status']=='SUFFICIENT' and temp['ci95_halfwidth'] is not None and abs(temp['mean']/spec['kBT_star']-1)+temp['ci95_halfwidth']/spec['kBT_star']<=tol['temperature_relative_error'] else 'INCONCLUSIVE_OR_FAILED',
                  all_sampled_finite=True,max_speed_star=float(np.max(moments['max_speed_star'])),
                  legacy_short_window_status='AVAILABLE_DURATION_ONLY' if physical_span>=c.get('legacy_windows',{}).get('short_window_s',0.0002441406727828746) else 'WINDOW_INSUFFICIENT',
                  legacy_long_window_status='AVAILABLE_DURATION_ONLY' if physical_span>=c.get('legacy_windows',{}).get('long_window_s',0.0004882813455657492) else 'WINDOW_INSUFFICIENT',
                  legacy_vessel_gates_executed=False,
                  series={'time_star':t.tolist(),'kBT_thermal_star':moments['kBT_thermal_star'].tolist(), 'kBT_raw_star':moments['kBT_raw_star'].tolist(),
                          'pressure_time_star':native['time'].tolist(),'pressure_star':p.tolist(), 'native_COM_corrected_kBT_star':knative.tolist()})
    if task['kind']=='flow':
        average=np.mean(by[sel],axis=0)
        fit=fit_viscosity(y,average,task['force_star'],spec['m_star'],Ly,c['space']['output_y_bin_star'])
        slope_stats=block_statistics(ts,slope_signal[sel],min_duration=min_block,min_blocks=tol['min_blocks'])
        b=slope_stats['block_samples']; n=slope_stats['block_count']
        fits=[fit_viscosity(y,np.mean(by[sel[i*b:(i+1)*b]],axis=0),task['force_star'],spec['m_star'],Ly,c['space']['output_y_bin_star']) for i in range(n)]
        nus=[z['nu_star'] for z in fits]; finite_nus=[z for z in nus if z is not None]
        # Invert slope CI (nonlinear reciprocal) rather than pretending Gaussian nu.
        h=slope_stats.get('ci95_halfwidth'); slope=slope_stats['mean']; ci=None
        if h is not None and slope-h>0:
            ci=[task['force_star']/(spec['m_star']*(slope+h)),task['force_star']/(spec['m_star']*(slope-h))]
        fit.update(block_nu_star=nus,block_times_star=slope_stats['block_times_star'],slope_statistics=slope_stats,ci95_star=ci,
                   measured_profile_star=average.tolist(),y_star=y.tolist(),block_valid_count=len(finite_nus))
        if fit['nu_star'] is not None:
            fit['nu_si']=units.to_si(fit['nu_star'],'kinematic_viscosity');fit['mu_si']=fit['nu_si']*units.to_si(result['actual_n_star']*spec['m_star'],'mass_density')
            fit['ci95_si']=[units.to_si(z,'kinematic_viscosity') for z in ci] if ci else None
            fit['relaxation_measured_star']=Ly**2/(4*math.pi**2*fit['nu_star'])
            # Reject if measured relaxation invalidates chosen warmup/block duration.
            fit['measured_relaxation_covered']=bool(start>=tol['warmup_relaxation_multiples']*fit['relaxation_measured_star'] and b*(ts[1]-ts[0])>=tol['block_relaxation_multiples']*fit['relaxation_measured_star'])
            fit['sampling_status']='SUFFICIENT' if slope_stats['status']=='SUFFICIENT' and ci and fit['measured_relaxation_covered'] else 'WINDOW_INSUFFICIENT'
        result['viscosity']=fit
    return result

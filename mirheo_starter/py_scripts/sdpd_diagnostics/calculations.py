"""CPU audits of recorded observables, with explicit sample and phase semantics."""
import math
from pathlib import Path
import numpy as np
from py_scripts.fluid_physics.analysis import (read_csv,periodic_shape,fit_viscosity,
    block_statistics,correlation_time_samples,t95,native_pressure)


from .observations import cuda_layout


def recover_profile_moments(counts,means,sum_v2,mass):
    counts=np.asarray(counts,float);means=np.asarray(means,float);sum_v2=np.asarray(sum_v2,float)
    if np.any(counts<0) or means.shape!=(*counts.shape,3):raise ValueError('INVALID_BIN_MOMENTS')
    live=counts>0;N=np.sum(counts,axis=-1);B=np.sum(live,axis=-1)
    if np.any(N<=1) or np.any(N<=B):raise ValueError('INSUFFICIENT_THERMAL_DOF')
    means=np.where(live[...,None],means,0.)
    if not np.isfinite(means).all() or not np.isfinite(sum_v2).all():raise ValueError('NONFINITE_BIN_MOMENTS')
    mean=np.sum(counts[...,None]*means,axis=-2)/N[...,None];s2=np.sum(sum_v2,axis=-1)
    raw=mass*s2/(3*N);com=mass*(s2-N*np.sum(mean**2,axis=-1))/(3*(N-1))
    local=mass*(s2-np.sum(counts*np.sum(means**2,axis=-1),axis=-1))/(3*(N-B))
    return {'raw':raw,'COM':com,'coarse':local,'mean':mean,'N':N,'occupied_bins':B}


def mass_field_semantics(spec,worker_text,moments):
    compact=''.join(worker_text.split())
    if "moments.writerow([done,done*dt,N,N*spec['m_star']" in compact:
        meaning='SYSTEM_TOTAL_MASS';expected=moments['N']*spec['m_star']
    elif 'moments.writerow([done,done*dt,N,mass,' in compact:
        meaning='PARTICLE_MASS';expected=np.full(len(moments),spec['m_star'])
    else:raise ValueError('UNVERIFIED_WORKER_MASS_FIELD_VERSION')
    return {'meaning':meaning,'matches_executed_formula':bool(np.allclose(moments['mass_star'],expected,rtol=0,atol=0)),
            'particle_mass_star':spec['m_star'],'system_mass_star':float(moments['N'][0]*spec['m_star'])}


def profile_arrays(directory,spec):
    m=read_csv(Path(directory)/'moments.csv');p=read_csv(Path(directory)/'profile_samples.csv')
    bins=int(round(spec['domain_star'][1]/spec['y_bin_star']))
    if len(p)!=len(m)*bins:raise ValueError('PROFILE_SAMPLE_COUNT_MISMATCH')
    time=p['time_star'].reshape(-1,bins);step=p['step'].reshape(-1,bins)
    if not np.allclose(time,m['time_star'][:,None],rtol=0,atol=1e-12) or not np.all(step==m['step'][:,None]):
        raise ValueError('PROFILE_PHASE_MISMATCH')
    if not np.all(p['bin'].reshape(-1,bins)==np.arange(bins)):raise ValueError('PROFILE_BIN_ORDER_MISMATCH')
    return m,p,p['count'].reshape(-1,bins),np.stack([p[n].reshape(-1,bins) for n in ['ux','uy','uz']],axis=-1),p['sum_v2'].reshape(-1,bins)


def window_audit(t,thermal,signal,tau_prior,tol,fractions):
    """Return every fixed candidate; reproduce the existing earliest-cut rule."""
    rows=[];chosen=None
    for frac in fractions:
        start=max(t[-1]*frac,tol['warmup_relaxation_multiples']*tau_prior)
        sel=np.flatnonzero(t>=start)
        if len(sel)<16:continue
        half=len(sel)//2;a,b=sel[:half],sel[half:]
        td=abs(float(np.mean(thermal[b])-np.mean(thermal[a])))/max(abs(float(np.mean(thermal[b]))),1e-30)
        sd=abs(float(np.mean(signal[b])-np.mean(signal[a])))/max(abs(float(np.mean(signal[b]))),1e-30)
        rec={'cut_fraction':frac,'start_star':float(t[sel[0]]),'end_star':float(t[-1]),'sample_count':len(sel),
             'temperature_mean_all':float(np.mean(thermal[sel])),'temperature_first_half':float(np.mean(thermal[a])),
             'temperature_second_half':float(np.mean(thermal[b])),'thermal_drift':td,'signal_drift':sd,
             'legacy_stationarity_pass':max(td,sd)<=tol['max_stationarity_drift'],
             'target_temperature_margin':tol['temperature_relative_error'],'temperature_ACF_samples':correlation_time_samples(thermal[sel]),
             'ACF_caution':'A drifting record is not a stationary-process correlation estimate.'}
        if frac in (.2,.3,.4,.5) and chosen is None and rec['legacy_stationarity_pass']:chosen=len(rows)
        rows.append(rec)
    if chosen is None:
        legacy=[i for i,r in enumerate(rows) if r['cut_fraction'] in (.2,.3,.4,.5)]
        chosen=legacy[-1] if legacy else None
    for i,r in enumerate(rows):r['selected_by_original_rule']=i==chosen
    return rows


def sampling_row(task,spec,kind,stats,tau_prior,tol,u,wall_s):
    minimum=tol['block_relaxation_multiples']*tau_prior
    spacing=stats.get('sample_interval_star');b=stats.get('block_samples');n=stats.get('block_count',0)
    return {'task_id':task,'observable':kind,'desired_steps':spec.get('desired_steps',spec['steps']),
        'allocated_steps':spec['steps'],'dt_star':spec['dt_star'],'planned_duration_star':spec['steps']*spec['dt_star'],
        'selected_start_star':stats.get('interval_star',[None,None])[0],
        'selected_end_star':stats.get('interval_star',[None,None])[1],
        'selected_samples':stats.get('sample_count',0),'sample_interval_star':spacing,
        'sample_interval_s':spacing*u.t0 if spacing else None,'tau_samples':stats.get('correlation_time_samples'),
        'tau_star':stats['correlation_time_samples']*spacing if spacing else None,
        'minimum_block_star':minimum,'minimum_block_samples':math.ceil(minimum/spacing) if spacing else None,
        'ACF_block_samples':math.ceil(5*stats['correlation_time_samples']) if spacing else None,
        'final_block_samples':b,'block_duration_star':stats.get('block_duration_star'),
        'block_duration_s':stats['block_duration_star']*u.t0 if b else None,
        'discarded_tail_samples':stats.get('discarded_tail_samples'),'block_count':n,
        'mean_before_tail_alignment':stats.get('all_sample_mean',stats.get('mean')),
        'mean_complete_blocks':float(np.mean(stats['block_means'])) if n else None,
        'ci95_halfwidth':stats.get('ci95_halfwidth'),'status':stats['status'],'task_wall_s':wall_s,
        'warning':'Selected interval is descriptive; a drift-based ACF is not proof of stationary independent samples.'}


def temperature_audit(directory,spec,old,corrected,u,fractions,tol):
    m,p,counts,means,v2=profile_arrays(directory,spec);t=m['time_star'];N=int(m['N'][0]);mass=spec['m_star']
    recovered=recover_profile_moments(counts,means,v2,mass)
    native=read_csv(Path(directory)/'native_stats.csv');virial=read_csv(Path(directory)/'pressure/pv.csv')
    pressure,knative=native_pressure(native,virial,mass,float(np.prod(spec['domain_star'])),spec['dt_star'])
    # Native Stats samples steps 1, every+1, ...; the worker samples every, 2*every, ... .
    native_step=np.rint(native['time']/spec['dt_star']).astype(int)
    native_expected=np.arange(len(native))*spec['native_stats_every']+1
    # Late six-digit CSV timestamps may hide one half-step; source schedule is authoritative.
    quantum=np.where(abs(native['time'])>0,10.**(np.floor(np.log10(np.maximum(abs(native['time']),1e-300)))-5),0.)
    schedule_ok=bool(np.all(abs(native['time']-native_expected*spec['dt_star'])<=quantum*.5+1e-10))
    index={int(s):i for i,s in enumerate(m['step'])};pairs=[(j,index[s-1]) for j,s in enumerate(native_expected) if s-1 in index]
    delta=np.array([knative[j]-m['kBT_COM_star'][i] for j,i in pairs])
    differences={'raw_max_abs':float(np.max(abs(recovered['raw']-m['kBT_raw_star']))),
        'COM_max_abs':float(np.max(abs(recovered['COM']-m['kBT_COM_star']))),
        'coarse_max_abs':float(np.max(abs(recovered['coarse']-m['kBT_bin_coarse_star']))) if 'kBT_bin_coarse_star' in m.dtype.names else None,
        'native_minus_worker_adjacent_step_mean':float(np.mean(delta)) if len(delta) else None,
        'native_minus_worker_adjacent_step_RMS':float(np.sqrt(np.mean(delta**2))) if len(delta) else None,
        'native_worker_are_same_instant':False,'native_worker_offset_steps':1,
        'source_independence':'CPU reconstruction uses stored bin moments from the same velocity getters, not independently saved old particle velocities.'}
    thermal=m['kBT_thermal_star'];Ly=spec['domain_star'][1]
    shape=periodic_shape(p['y_star'][:counts.shape[1]],Ly,spec['y_bin_star'])
    signal=means[:,:,0]@shape/np.dot(shape,shape) if spec['task']['kind']=='flow' else thermal
    tau=Ly**2/(4*math.pi**2*spec['candidate']['nu_prior_star'])
    windows=window_audit(t,thermal,signal,tau,tol,fractions)
    coarse=recovered['coarse'];row={'task_id':spec['task']['id'],'N':N,'mass_semantics':mass_field_semantics(spec,(Path(directory)/'gpu_worker.py').read_text(),m),
        'DOF_relative_correction':1/(N-1),'initial_native_kBT_COM':float(knative[0]),
        'first_worker_kBT':float(thermal[0]),'peak_kBT':float(np.max(thermal)),'peak_time_star':float(t[np.argmax(thermal)]),
        'last_sample_kBT':float(thermal[-1]),'last_10pct_mean':float(np.mean(thermal[t>=.9*t[-1]])),
        'last_10pct_range':[float(np.min(thermal[t>=.9*t[-1]])),float(np.max(thermal[t>=.9*t[-1]]))],
        'original_selected_mean':old.get('temperature',{}).get('mean'),'corrected_selected_mean':corrected.get('temperature',{}).get('mean'),
        'original_temperature_status':old.get('temperature_status'),'max_COM_speed':float(np.max(np.linalg.norm(recovered['mean'],axis=-1))),
        'coarse_vs_COM_max_abs':float(np.max(abs(coarse-recovered['COM']))),'empty_bin_records':int(np.count_nonzero(counts==0)),
        'min_bin_count':int(counts.min()),'native_schedule_verified':schedule_ok,'independent_moment_recovery':differences,
        'velocity_component_temperatures':None,'component_temperature_reason':'Old moments store total sum_v2, not sum_vx2/sum_vy2/sum_vz2; old particle velocity snapshots were not saved.',
        'old_particle_outlier_energy_distribution':None,'old_particle_outlier_reason':'Only max_speed and scalar/bin moments exist; the full velocity distribution cannot be recovered.',
        'windows':windows,'series':{'time_star':t.tolist(),'thermal':thermal.tolist(),'COM':m['kBT_COM_star'].tolist(),'coarse':coarse.tolist(),
        'native_time_star':native['time'].tolist(),'native_COM':knative.tolist(),'paired_time_star':[float(t[i]) for j,i in pairs],
        'adjacent_step_native_minus_worker':delta.tolist(),'pressure_time_star':native['time'].tolist(),'mechanical_pressure_pa':(pressure*u.scales['pressure']).tolist(),
        'spatial_coarse_temperature':(mass*(v2-counts*np.sum(means**2,axis=-1))/(3*np.maximum(counts-1,1))).tolist()}}
    if (Path(directory)/'density_samples.csv').exists():
        d=read_csv(Path(directory)/'density_samples.csv');row['series'].update(density_time_star=d['density_time_star'].tolist(),
            kernel_mean=d['kernel_number_mean'].tolist(),kernel_cv=(d['kernel_number_std']/d['kernel_number_mean']).tolist())
    return row


def profile_audit(directory,spec,old,corrected,u):
    if spec['task']['kind']!='flow':return None
    m,p,counts,means,v2=profile_arrays(directory,spec);nb=counts.shape[1];y=p['y_star'][:nb]
    start=old['statistics_interval_star'][0];use=m['time_star']>=start;v=means[use,:,0];weights=counts[use]
    equal,weighted=average_profile(v,weights)
    f=fit_viscosity(y,equal,spec['task']['force_star'],spec['m_star'],spec['domain_star'][1],spec['y_bin_star'])
    residual=equal-np.array(f['fit_profile_star']);h=spec['domain_star'][1]/2
    distance=np.minimum(np.mod(y,h),h-np.mod(y,h));near=distance<=spec['rc_star']
    return {'task_id':spec['task']['id'],'y_star':y.tolist(),'measured_all_window_star':equal.tolist(),
      'target_star':old['viscosity']['target_profile_star'],'old_fit_star':f['fit_profile_star'],
      'corrected_fit_star':corrected['viscosity']['fit_profile_star'],'corrected_measured_star':corrected['viscosity']['measured_profile_star'],
      'residual_star':residual.tolist(),'particle_weighted_profile_star':weighted.tolist(),
      'weighted_vs_equal_RMS_star':float(np.sqrt(np.mean((weighted-equal)**2))),
      'near_force_planes_RMS_star':float(np.sqrt(np.mean(residual[near]**2))),
      'interior_RMS_star':float(np.sqrt(np.mean(residual[~near]**2))),
      'near_force_planes_residual_energy_fraction':float(np.sum(residual[near]**2)/np.sum(residual**2)),
      'rc_over_half_channel':spec['rc_star']/h,'rc_over_bin_width':spec['rc_star']/spec['y_bin_star'],
      'bin_count_min':int(weights.min()),'empty_bin_records':int(np.sum(weights==0)),
      'mean_bin_density_relative_to_global':(weights.mean(axis=0)/(float(m['N'][0])/nb)).tolist(),
      'old_nu_si':old['viscosity']['nu_si'],'corrected_nu_si':corrected['viscosity']['nu_si'],
      'old_RMS_relative':old['viscosity']['relative_rms'],'corrected_RMS_relative':corrected['viscosity']['relative_rms'],
      'valid_temporal_blocks':corrected['viscosity']['block_valid_count'],
      'noise_explains_residual':None,'noise_reason':'No trustworthy stationary multi-block profile CI in this record; spatial residuals are descriptive, not proof of continuum/discrete model bias.'}


def shared_reference_covariance(variances):
    """Cov(P_low-P_ref, P_high-P_ref) for independent initialized runs."""
    low,ref,high=map(float,variances)
    if min(low,ref,high)<0:raise ValueError('NEGATIVE_VARIANCE')
    return np.array([[low+ref,ref],[ref,high+ref]])


def matched_physical_interval(first,second,dt_first,dt_second):
    t1=np.asarray(first)*dt_first;t2=np.asarray(second)*dt_second
    lower=max(float(t1[0]),float(t2[0]));upper=min(float(t1[-1]),float(t2[-1]))
    if upper<=lower:raise ValueError('NO_MATCHED_PHYSICAL_INTERVAL')
    return lower,upper


def matched_flow_audit(base,other,u):
    """Compare both recorded flows on the same physical sample windows."""
    records=[]
    for row in (base,other):
        spec=row['parameters'];m,p,n,v,_=profile_arrays(row['directory'],spec)
        records.append((spec,m,p,n,v))
    end=min(r[1]['time_star'][-1] for r in records);windows=[]
    for lo,hi in [(0,end/2),(end/2,end)]:
        pair=[]
        for spec,m,p,n,v in records:
            mask=(m['time_star']>lo)&(m['time_star']<=hi+1e-12)
            velocity=average_profile(v[mask,:,0],n[mask])[0]
            fit=fit_viscosity(p['y_star'][:n.shape[1]],velocity,spec['task']['force_star'],spec['m_star'],spec['domain_star'][1],spec['y_bin_star'])
            pair.append({'task_id':spec['task']['id'],'samples':int(mask.sum()),'temperature_mean':float(m['kBT_thermal_star'][mask].mean()),
              'nu_fit_si':float(u.to_si(fit['nu_star'],'kinematic_viscosity')),'relative_RMS':fit['relative_rms'],
              'force_star':spec['task']['force_star'],'measured_profile_star':velocity.tolist(),
              'velocity_per_particle_force':(velocity/spec['task']['force_star']).tolist()})
        windows.append({'interval_star':[lo,float(hi)],'records':pair})
    return {'task_id':other['task_id'],'joint_interval_star':[0,float(end)],'joint_duration_s':float(end*u.t0),
      'windows':windows,'status':'TRANSIENT_COMPARISON_ONLY','other_execution':other['execution']['status'],
      'reason':'Same physical stages; both records retain startup relaxation. Force-normalized fits are descriptive, not steady convergence or independent replicates.'}


def average_profile(velocity,counts):
    v=np.asarray(velocity,float);n=np.asarray(counts,float);live=n>0
    if v.shape!=n.shape or np.any(live.sum(axis=0)==0):raise ValueError('PROFILE_BIN_WITHOUT_OBSERVATIONS')
    clean=np.where(live,v,0.)
    if not np.isfinite(clean).all():raise ValueError('NONFINITE_LIVE_PROFILE')
    return clean.sum(axis=0)/live.sum(axis=0),(clean*n).sum(axis=0)/n.sum(axis=0)


def ou_discrete_temperature_ratio(damping_rate,dt):
    x=float(damping_rate)*float(dt)
    if not 0<=x<2:raise ValueError('EXPLICIT_OU_UNSTABLE')
    return 1/(1-x/2)

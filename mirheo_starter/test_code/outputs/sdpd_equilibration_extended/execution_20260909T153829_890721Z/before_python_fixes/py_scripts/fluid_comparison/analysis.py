"""Reanalyse immutable native runs with explicit method/candidate/group identities."""
from pathlib import Path
import time
import numpy as np
from py_scripts.fluid_physics.common import read_json,sha256_file
from py_scripts.fluid_physics.analysis import analyze_task,read_csv,periodic_shape,native_pressure
from .models import snapshot_density


def checked_run(directory):
    d=Path(directory)
    for name,h in read_json(d/'output_sha256.json').items():
        p=(d/name).resolve()
        if not p.is_relative_to(d.resolve()) or sha256_file(p)!=h:raise ValueError('RAW_RUN_HASH_MISMATCH '+str(p))
    return d


def analyze_run(directory,c,u,case,*,historical=False):
    d=checked_run(directory);r=analyze_task(d,c,u);spec=read_json(d/'actual_parameters.json')
    r['source_category']='HISTORICAL_REFERENCE' if historical else 'MEASURED'
    r['execution']=read_json(d/'execution.json')
    if 'ended_at_utc' not in r['execution']:
        r['execution']['ended_at_utc']=r['execution']['started_at_utc']
        r['execution']['started_at_utc']=read_json(d/'process_identity.json')['started_at']
        r['execution']['annotation']='Historical UTC field correction only; original files and monotonic duration preserved.'
    try:m=read_csv(d/'moments.csv')
    except ValueError as exc:r['read_error']=str(exc);return r
    r['actual_N']=int(m['N'][0]);r['actual_n_star']=float(m['N'][0]/np.prod(spec['domain_star']))
    r['actual_rho_star']=r['actual_n_star']*spec['m_star'];r['density_si']=u.to_si(r['actual_rho_star'],'mass_density')
    r['mass_conserved']=bool(np.all(m['N']==m['N'][0]))
    r.setdefault('series',{}).update(time_star=m['time_star'].tolist(),kBT_thermal_star=m['kBT_thermal_star'].tolist(),
                                   kBT_raw_star=m['kBT_raw_star'].tolist(),kBT_COM_star=m['kBT_COM_star'].tolist())
    try:
        native=read_csv(d/'native_stats.csv');virial=read_csv(d/'pressure/pv.csv')
        pressure,_=native_pressure(native,virial,spec['m_star'],float(np.prod(spec['domain_star'])),spec['dt_star'])
        r['pressure_channels_verified']=True
        r['series'].update(pressure_time_star=native['time'].tolist(),pressure_star=pressure.tolist())
        r['pressure_time_audit']={'expected_offset_star':spec['dt_star'],'observed_offset_range_star':[float(np.min(native['time']-virial['time'])),float(np.max(native['time']-virial['time']))],
          'tolerance':'Half a 6-significant-digit CSV quantum for each timestamp plus 1e-10; late labels may not individually resolve a half-micro-step. Early resolved labels and pinned callbacks constrain the offset.'}
    except (ValueError,OSError) as exc:r['pressure_channels_verified']=False;r['pressure_time_audit']={'error':str(exc)}
    start=r.get('statistics_interval_star',[None])[0]
    sel=m['time_star']>=start if start is not None else np.zeros(len(m),dtype=bool)
    r['initialization_temperature_range_star']=[float(np.min(m['kBT_thermal_star'])),float(np.max(m['kBT_thermal_star']))]
    if 'kBT_bin_coarse_star' in m.dtype.names:
        r['series']['kBT_bin_coarse_star']=m['kBT_bin_coarse_star'].tolist()
        if np.any(sel):
            difference=float(abs(np.mean(m['kBT_bin_coarse_star'][sel])-np.mean(m['kBT_thermal_star'][sel]))/spec['kBT_star'])
            r['temperature_refinement']={'coarse_bin_star':spec['y_bin_star'],'fine_bin_star':spec['y_bin_star']/2,
               'relative_difference_to_target':difference,'status':'PASS_PROPOSED' if difference<=c['proposed_additional']['temperature_bin_refinement_difference'] else 'FAIL_PROPOSED',
               'scope':'Flow: independently fitted instantaneous bin means, with each fitted mean accounted in degrees of freedom. Equilibrium uses COM temperature.'}
    else:r['temperature_refinement']={'status':'HISTORICAL_COARSE_ONLY','scope':'No stored finer bins/coordinates; do not fabricate a refined measurement.'}
    v=r.get('viscosity')
    if v and v.get('nu_star'):
        target_star=u.to_star(case['kinematic_viscosity_m2_s'],'kinematic_viscosity')
        v['target_profile_star']=(spec['task']['force_star']/(spec['m_star']*target_star)*periodic_shape(v['y_star'],spec['domain_star'][1],spec['y_bin_star'])).tolist()
        v['input_mu_star']=spec['candidate'].get('viscosity_mu_star')
        v['measured_mu_star']=v['nu_star']*r['actual_rho_star'];v['source_category']='FITTED'
        v['target_relative_error']=abs(v['nu_si']/case['kinematic_viscosity_m2_s']-1)
        # Estimate unresolved smooth shear using the measured fit, never subtract an arbitrary target field.
        y=np.array(v['y_star']);h=spec['domain_star'][1]/2;z=np.where(y<h,y,y-h)
        shear=spec['task']['force_star']/(spec['m_star']*v['nu_star'])*(h-2*z)/2
        v['estimated_coarse_bin_shear_temperature_residue_star']=float(spec['m_star']*np.mean(shear**2)*spec['y_bin_star']**2/36)
    complete=read_json(d/'worker_completion.json') if (d/'worker_completion.json').exists() else {}
    elapsed=r['execution']['elapsed_monotonic_s'];steps=int(complete.get('steps',m['step'][-1]));duration=steps*spec['dt_star']*u.t0
    compute=complete.get('compute_wall_s');sample=complete.get('snapshot_wall_s')
    perf={'total_wall_s':elapsed,'actual_steps':steps,'actual_physical_duration_s':duration,
          'statistics_physical_duration_s':r.get('statistics_duration_s'),'total_wall_s_per_step':elapsed/steps,
          'wall_s_per_physical_s':elapsed/duration,'compute_wall_s':compute,'sampling_transfer_output_s':sample,
          'worker_setup_wall_s':complete.get('setup_wall_s'),'first_chunk_compute_s':complete.get('first_chunk_compute_s'),
          'synchronized_steady_compute_s_per_step':complete.get('steady_compute_wall_s',0)/complete['steady_steps'] if complete.get('steady_steps',0)>0 else None,
          'timing_scope':'synchronized native run blocks include scheduler, MPI and native plugins; sampling includes CPU transfer and files' if not historical else 'Historical compute timer lacks separate explicit completion event; combined compute+sampling is bounded by synchronized getters. Compare only total cost as historical reference.',
          'precision':spec['precision'],'N':r['actual_N'],'domain_star':spec['domain_star'],'rc_star':spec['rc_star'],
          'density_kernel':'WendlandC2' if complete.get('density_enabled') else None,'sample_interval_physical_s':spec['snapshot_every']*spec['dt_star']*u.t0,
          'device_sampled_peak_used_MiB':r['execution'].get('device_sampled_peak_used_MiB'),
          'memory_scope':r['execution'].get('memory_scope'),'time_to_qualified_solution_s':None}
    mem=[s['host_memory_bytes']['MemAvailable']/1024**2 for s in r['execution'].get('resource_samples',[]) if s.get('host_memory_bytes')]
    perf['host_min_available_MiB']=min(mem) if mem else None
    if compute is not None and sample is not None:
        perf['outside_worker_setup_compute_sampling_s']=elapsed-compute-sample-complete.get('setup_wall_s',0)
    if (d/'chunk_timing.csv').exists() and start is not None:
        chunks=read_csv(d/'chunk_timing.csv');warm=chunks['time_star']<start
        perf['equilibration_blocks_wall_s']=float(np.sum(chunks['compute_sync_s'][warm]+chunks['sampling_transfer_output_s'][warm]))
        perf['statistics_blocks_compute_s']=float(np.sum(chunks['compute_sync_s'][~warm]))
        perf['statistics_blocks_sampling_s']=float(np.sum(chunks['sampling_transfer_output_s'][~warm]))
    r['performance']=perf
    if spec['candidate'].get('method')=='SDPD':
        density=read_csv(d/'density_samples.csv');r['series']['density_time_star']=density['density_time_star'].tolist()
        for name in ['kernel_number_mean','kernel_number_std','kernel_mass_mean','EOS_model_mean_star']:
            r['series'][name]=density[name].tolist()
        use=density['density_time_star']>=start if start is not None else np.zeros(len(density),dtype=bool)
        if np.count_nonzero(use)>=8:
            means=density['kernel_number_mean'][use];std=density['kernel_number_std'][use];mid=len(means)//2
            cv=float(np.mean(std/means));drift=float(abs(np.mean(means[mid:])/np.mean(means[:mid])-1))
            r['kernel_density']={'global_number_density':r['actual_n_star'],'global_mass_density':r['actual_rho_star'],
              'local_number_mean':float(means.mean()),'local_mass_mean':float(means.mean()*spec['m_star']),
              'local_number_bias_relative_to_global':float(means.mean()/r['actual_n_star']-1),
              'local_cv_mean':cv,'mean_drift_relative':drift,
              'status':'PASS_PROPOSED' if cv<=c['proposed_additional']['kernel_density_cv_max'] and drift<=c['proposed_additional']['kernel_density_time_drift_max'] else 'FAIL_PROPOSED',
              'model_EOS_pressure_mean_star':float(np.mean(density['EOS_model_mean_star'][use])),
              'pressure_category':'DERIVED_FROM_MEASURED_LOCAL_KERNEL_DENSITY; not independent mechanical pressure'}
        else:r['kernel_density']={'status':'WINDOW_INSUFFICIENT'}
        checks=[];audit_started=time.monotonic()
        paths=sorted(d.glob('snapshot_*_positions.npy'))
        for p in paths[-2:]:
            meta=read_json(p.with_name(p.name.replace('_positions.npy','.json')))
            observed=np.load(p.with_name(p.name.replace('_positions.npy','_kernel_number.npy')),allow_pickle=False)
            calculated=snapshot_density(np.load(p,allow_pickle=False),spec['domain_star'],spec['rc_star'])
            error=float(np.max(abs(calculated-observed)/np.maximum(abs(observed),1e-30)))
            checks.append({'snapshot':p.name,'phase':meta,'max_relative_error':error,
                           'status':'PASS' if error<=c['proposed_additional']['kernel_CPU_relative_check'] else 'FAIL'})
        r['kernel_CPU_audit']={'checks':checks,'cpu_wall_s':time.monotonic()-audit_started,'category':'DERIVED_INDEPENDENT_SNAPSHOT_CHECK',
                               'status':'PASS' if checks and all(z['status']=='PASS' for z in checks) else 'UNVERIFIED_OR_FAILED'}
    return r

"""Analyze actual owned-fluid samples; do not fill missing measurements."""
from pathlib import Path
import math
import numpy as np
from py_scripts.fluid_physics.common import read_json,sha256_file,PROJECT_ROOT
from py_scripts.fluid_physics.analysis import t95
from .physics import target,profile_metrics,valid_blocks,qualified_ratio
from .workflow import verify_run,paths


def read_profiles(path,bins):
    a=np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True))
    if len(a)%bins or not all(np.isfinite(a[n]).all() for n in a.dtype.names):raise ValueError('NONFINITE_OR_INCOMPLETE_PROFILE')
    for offset in range(0,len(a),bins):
        row=a[offset:offset+bins]
        if len(np.unique(row['step']))!=1 or not np.array_equal(row['bin'],np.arange(bins)):raise ValueError('DUPLICATED_BINS_OR_HALOS')
    return a,a['time_si'][::bins],a['ux_si'].reshape(-1,bins),a['count'].reshape(-1,bins)


def timing_metrics(completion,execution,timings,p,c,task=None):
    n=completion['actual_steps'];elapsed=execution['elapsed_monotonic_s'];duration=completion['actual_time_si']
    rows=np.atleast_1d(timings);core=float(np.sum(rows['compute_max_rank_s']))
    sampling=float(np.sum(rows['sampling_max_rank_s']));rates=rows['compute_max_rank_s']/rows['chunk_steps']
    phases={}
    warm_end=task['warmup_steps']*task['dt_si'] if task else c['sampling']['unforced_preparation_si']
    for name,lo,hi in [('preparation',0,warm_end),('flow_establishment',warm_end,c['sampling']['flow_establishment_end_si']),('formal_statistics',c['sampling']['flow_establishment_end_si'],float('inf'))]:
        mask=(rows['time_si']>lo)&(rows['time_si']<=hi)
        phases[name]={'steps':int(rows['chunk_steps'][mask].sum()),'compute_s':float(rows['compute_max_rank_s'][mask].sum()),'sampling_s':float(rows['sampling_max_rank_s'][mask].sum())}
    return dict(total_wall_s=elapsed,setup_s=completion.get('setup_s'),geometry_s=completion.get('geometry_s'),
        compute_s=core,sampling_transfer_reduction_output_s=sampling,
        launch_setup_exit_and_unclassified_s=elapsed-core-sampling,
        compute_ms_per_step=core/n*1000,workflow_ms_per_step=elapsed/n*1000,
        core_wall_s_per_us=core/duration*1e-6,workflow_wall_s_per_us=elapsed/duration*1e-6,
        physical_us_per_wall_s=duration*1e6/elapsed,
        chunk_ms_per_step_median=float(np.median(rates)*1000),chunk_ms_per_step_p10=float(np.quantile(rates,.1)*1000),chunk_ms_per_step_p90=float(np.quantile(rates,.9)*1000),
        compute_timing_scope='Synchronized native steps; MPI maximum per chunk. No build or download costs.',
        sampling_timing_scope='Host transfer (Mirheo), spatial moments, MPI reduction and CSV; timing CSV write and orchestration remain in total wall time.',
        phases=phases,time_to_qualified_solution_s=None)


def analyze_run(d,c,p):
    d=Path(d);verify_run(d);task=read_json(d/'task.json');ex=read_json(d/'execution.json')
    backend=task['backend'];is_lbm=backend=='HemoCell';role=task.get('role','main')
    base=dict(backend=backend,method='HemoCell/Palabos纯流体路径：D3Q19 Guo BGK' if is_lbm else 'Mirheo SDPD / Wendland C2 / Linear EOS',
        case_id=c['case_id'],candidate_id=('lbm_N'+str(task['N']) if is_lbm else p['candidate']['id']),
        task_id=task.get('task_id',d.name),role=role,repetition=task.get('repetition',0),source_commit=task['source_commit'],
        binary_hash=task['binary_hash'],precision='float64' if is_lbm else 'float32',rank_count=task.get('rank_count',2),
        hardware='i7-13700HX CPU / WSL2' if is_lbm else 'RTX 4060 Laptop GPU / 1 compute + 1 postprocess rank',
        source_category='MEASURED',directory=str(d),physical_parameters=p,dt_si=task['dt_si'],
        resolution={'dx_si':task.get('dx_si'),'nodes':task.get('lattice_nodes'),'particle_spacing_si':None if is_lbm else p['mapping']['particle_spacing_m'],
                    'kernel_radius_si':None if is_lbm else p['mapping']['interaction_cutoff_m'],'bins':p['bins'],'bin_width_si':p['bin_width_si']},
        qualified=False,time_to_qualified_solution_s=None,budget_status=ex['status'],limitations=[])
    if not (d/'completion.json').exists():
        console=(d/'console.log').read_text();preimport="No module named 'mpi4py'" in console
        host=read_json(d/'host_memory.json')
        if not host.get('ranks'):
            host={**host,'coverage':'LAUNCHER_ONLY_RANK_RSS_NOT_MEASURED','sampled_peak_launcher_rss_bytes':host.get('sampled_peak_tree_rss_bytes'),
                  'sampled_peak_tree_rss_bytes':None,'scope':'Only the mpirun launcher was captured; no solver tree RSS measurement.'}
        return {**base,'actual_steps':0 if preimport else None,'actual_time_si':0 if preimport else None,'accuracy':None,'timing':{'total_wall_s':ex['elapsed_monotonic_s']},
                'sampling_status':'STARTUP_FAILED_BEFORE_SOLVER_IMPORT' if preimport else 'NOT_AVAILABLE','status':'FAILED_OR_STOPPED',
                'memory':{'host':host,'device_sampled_peak_used_MiB':ex.get('device_sampled_peak_used_MiB'),'device_scope':'Startup-only device-wide background sample; not a fluid solver memory measurement'},
                'limitations':["Missing mpi4py in the preserved virtual environment; 0 solver steps. Corrected to installed OpenMPI C ABI; explicit retry required." if preimport else 'Worker completion missing; inspect console.log. No automatic retry.']}
    comp=read_json(d/'completion.json');a,t,profiles,counts=read_profiles(d/'profiles.csv',p['bins'])
    tm=np.genfromtxt(d/'timings.csv',delimiter=',',names=True)
    base.update(actual_steps=comp['actual_steps'],actual_time_si=comp['actual_time_si'],completion=comp,
                timing=timing_metrics(comp,ex,tm,p,c,task),status='COMPLETED' if comp['completed'] else 'PARTIAL',
                actual_unforced_preparation_end_si=task['warmup_steps']*task['dt_si'])
    mem=read_json(d/'host_memory.json') if (d/'host_memory.json').exists() else {}
    if not mem.get('ranks'):
        mem={**mem,'coverage':'LAUNCHER_ONLY_RANK_RSS_NOT_MEASURED',
             'raw_reported_tree_peak_bytes':mem.get('sampled_peak_tree_rss_bytes'),
             'sampled_peak_launcher_rss_bytes':mem.get('sampled_peak_tree_rss_bytes'),'sampled_peak_tree_rss_bytes':None,
             'scope':'Original process-group sampler captured launcher only because OpenMPI rank groups differ. Raw record preserved; solver tree/rank RSS is unknown.'}
    base['memory']={'host':mem,'device_sampled_peak_used_MiB':ex.get('device_sampled_peak_used_MiB') if not is_lbm else None,
                    'device_scope':'Device-wide sampled peak includes desktop/other programs; never RSS' if not is_lbm else 'NOT_APPLICABLE_CPU_BACKEND'}
    y,truth=target(p);norm=np.linalg.norm(truth);half=p['bins']//2
    errs=np.linalg.norm(profiles-truth,axis=1)/norm
    base['series']={'time_si':t.tolist(),'upper_mean_si':profiles[:,half:].mean(axis=1).tolist(),
                    'lower_mean_si':profiles[:,:half].mean(axis=1).tolist(),'profile_relative_l2':errs.tolist(),
                    'temperature_K':None,'temperature_time_si':None}
    expected_count=task['lattice_nodes'] if is_lbm else p['mapping']['n_star']*np.prod(p['domain_star'])
    if not np.all(counts.sum(axis=1)==expected_count):raise ValueError('OWNED_FLUID_COUNT_MISMATCH')
    if is_lbm:
        if not np.all(counts==expected_count/p['bins']):raise ValueError('HALO_OR_NONUNIFORM_BIN_COUNTS')
        for key in ['dt_si','tau','dx_si']:
            if not math.isclose(comp[key],task[key],rel_tol=1e-12):raise ValueError('ACTUAL_LATTICE_PARAMETER_MISMATCH')
        if not math.isclose((comp['tau']-.5)*comp['dx_si']**2/(3*comp['dt_si']),p['nu_si'],rel_tol=1e-12):raise ValueError('ACTUAL_VISCOSITY_MAPPING_MISMATCH')
        correction=float(np.max(np.abs(a['ux_si']-a['raw_ux_si']-a['half_force_si'])))
        if correction>1e-12:raise ValueError('GUO_VELOCITY_CORRECTION_MISMATCH')
        base['velocity_half_force_check_max_si']=correction
        base['Mach_measured_max']=float(np.max(a['max_speed_si'])*task['dt_si']/task['dx_si']*math.sqrt(3))
        base['lattice_max_velocity']=base['Mach_measured_max']/math.sqrt(3)
        base['temperature_status']='NOT_APPLICABLE_ISOTHERMAL_LBM'
    else:
        mom=np.atleast_1d(np.genfromtxt(d/'moments.csv',delimiter=',',names=True))
        if not np.array_equal(mom['time_si'],t):raise ValueError('PROFILE_TEMPERATURE_PHASE_MISMATCH')
        for key in mom.dtype.names:
            values=mom[key][1:] if key in ['kernel_rho_mean_si','kernel_rho_cv'] else mom[key]
            if not np.isfinite(values).all():raise ValueError('NONFINITE_MOMENTS '+key)
        base['series'].update(temperature_K=mom['temperature_K'].tolist(),temperature_time_si=t.tolist())
        base['Mach_measured_max']=float(np.max(mom['max_speed_si'])/(p['candidate']['sound_speed_star']*p['mapping']['si_per_star']['velocity']))
    selected=np.flatnonzero((t>c['sampling']['flow_establishment_end_si'])&(t<=c['sampling']['formal_end_si']+task['dt_si']/2))
    # Exclude only an incomplete temporal tail for regular block statistics; no spatial bins removed.
    tail=[]
    if len(selected)>3:
        spacing=float(np.median(np.diff(t[selected])))
        if abs(t[selected[-1]]-t[selected[-2]]-spacing)>.02*spacing:tail=[int(selected[-1])];selected=selected[:-1]
    if role!='main' or len(selected)<8:
        base.update(accuracy={'last_profile':profile_metrics(p,profiles[-1]),'qualified':False},sampling_status='SMOKE_ONLY' if role=='smoke' else 'EXECUTION_COST_ONLY',
                    limitations=['这是执行成本，不是获得正确稳态解的耗时；部署冒烟只验证零细胞、非零步进、输出回读。' if role=='smoke' else '低 I/O 辅助成本段，采样间隔不同；不用于合格解速度比。'])
        if role=='memory_audit':base.update(sampling_status='RESOURCE_AUDIT_ONLY',limitations=['Supplemental fresh run at the same main-case schedule for descendant/rank RSS capture. Main timing repetitions remain separate.'])
        return base
    window=profiles[selected];tw=t[selected];mid=len(selected)//2
    drift=float(np.linalg.norm(window[:mid].mean(axis=0)-window[mid:].mean(axis=0))/norm)
    mean=window.mean(axis=0);ci=None;sample_status='DETERMINISTIC_NO_THERMAL_CI';stat={};used=len(selected)
    if not is_lbm:
        stats=[valid_blocks(tw,window[:,i],p,c) for i in range(p['bins'])]
        block=max(x.get('block_samples',len(selected)+1) for x in stats);n=len(selected)//block
        used=n*block if n else len(selected)
        mean=window[:used].mean(axis=0)
        sufficient=n>=c['criteria']['min_blocks'] and drift<=c['criteria']['stationarity_relative_drift']
        if sufficient:
            means=window[:n*block].reshape(n,block,p['bins']).mean(axis=1)
            ci=t95(n-1)*means.std(axis=0,ddof=1)/math.sqrt(n)
        stat={'block_samples':block,'block_count':n,'estimator_sample_count':used,'discarded_tail_samples':len(selected)-n*block,
              'per_bin_correlation_time_samples':[x.get('correlation_time_samples') for x in stats],
              'CI_method':'Complete nonoverlapping block means, >=5 autocorrelation times and >=2 viscous relaxation times; >=8 blocks and stationarity required'}
        sample_status='SUFFICIENT' if sufficient else 'INSUFFICIENT_INDEPENDENT_BLOCKS_OR_NONSTATIONARY'
        temperature=mom['temperature_K'][selected];tstat=valid_blocks(tw,temperature,p,c)
        tdrift=abs(float(temperature[:mid].mean()-temperature[mid:].mean()))/p['temperature_K']
        if tdrift>c['criteria']['stationarity_relative_drift']:tstat['ci95_halfwidth']=None
        temppass=tstat.get('ci95_halfwidth') is not None and abs(tstat['mean']/p['temperature_K']-1)+tstat['ci95_halfwidth']/p['temperature_K']<=c['criteria']['mirheo_temperature_relative_error_including_ci']
        base['temperature_status']='PASS_PROPOSED' if temppass else 'FAILED_OR_INCONCLUSIVE'
        base['temperature_statistics']={**tstat,'drift_relative_to_target':tdrift,'measured_mean_K':tstat['mean'],'target_K':p['temperature_K'],
            'bin_refinement_difference_relative':float(np.mean(np.abs(mom['temperature_K'][selected]-mom['temperature_coarse_K'][selected]))/p['temperature_K'])}
        base['density_diagnostics']={'kernel_mass_density_mean_si':float(mom['kernel_rho_mean_si'][selected].mean()),
            'kernel_mass_density_mean_min_si':float(mom['kernel_rho_mean_si'][selected].min()),'kernel_mass_density_mean_max_si':float(mom['kernel_rho_mean_si'][selected].max()),
            'kernel_mass_density_cv_mean':float(mom['kernel_rho_cv'][selected].mean()),'scope':'Native self-inclusive kernel mass density at pre-integration phase; distinct from global N*m/V and Eulerian bin particle counts.'}
    accuracy=profile_metrics(p,mean)
    rhobar=a['rho_si'].reshape(-1,p['bins'])[selected]
    global_error=float(np.max(abs(rhobar.mean(axis=1)/p['rho_si']-1)))
    local_variation=float(np.max(abs(rhobar/p['rho_si']-1)))
    if is_lbm:local_variation=float(max(np.max(abs(a['rho_min_si']/p['rho_si']-1)),np.max(abs(a['rho_max_si']/p['rho_si']-1))))
    accuracy.update(stationarity_relative_drift=drift,global_density_relative_error=global_error,
        local_density_relative_variation=local_variation,profile_ci95_si=ci.tolist() if ci is not None else None,
        statistical_error_status=sample_status,discretization_error_status='REFINEMENT_CHECK_PENDING' if is_lbm else 'SDPD_RESOLUTION_NOT_SCANNED',
        physical_model_error_status='NOT_QUANTIFIED; mean-flow screen does not validate real liquid, pressure outlets, walls, cells or thermal fluctuation equivalence',
        actual_statistics_interval_si=[float(tw[0]),float(tw[-1])],incomplete_time_tail_excluded=tail,
        profile_estimator_interval_si=[float(tw[0]),float(tw[used-1])],profile_estimator_samples=used,
        full_window_mean_profile_si=window.mean(axis=0).tolist(),full_window_descriptive_metrics=profile_metrics(p,window.mean(axis=0)),
        density_variation_scope='Maximum instantaneous Eulerian-bin density deviation (LBM also checks individual nodes); particle thermal count fluctuations retained, not a thermal-equivalence test.',
        space_bins_excluded=0,block_statistics=stat)
    criteria=c['criteria'];gates={}
    for key in ['profile_relative_l2','half_flow_relative_error','apparent_nu_relative_error','stationarity_relative_drift','global_density_relative_error']:
        gates[key]=accuracy.get(key) is not None and accuracy[key]<=criteria[key]
    gates['local_density_variation']=local_variation<=criteria['local_density_relative_variation']
    if not is_lbm:
        gates['valid_profile_CI']=ci is not None and np.linalg.norm(ci)/norm<=criteria['profile_ci_relative_l2']
        gates['temperature']=base['temperature_status']=='PASS_PROPOSED'
    gates['complete_fixed_window']=comp['completed'] and ex['status']=='COMPLETED'
    accuracy['gates']={k:bool(v) for k,v in gates.items()};accuracy['local_mean_flow_pass']=all(gates.values())
    base.update(accuracy=accuracy,sampling_status=sample_status)
    base['limitations']=['本次未验收真实压力出口及完整物性；0.263 Pa 历史门槛不是本周期速度初筛前置条件。',
        '相同平均间距不等于相同精度；LBM 格距与 SDPD 核半径、粒子间距分别记录。']
    return base


def analyze_all(c,frozen):
    _,runs=paths(c);rows=[]
    for backend in ['cpu','gpu']:
        ledger=runs/backend/'budget_ledger.json'
        if ledger.exists():
            for attempt in read_json(ledger)['attempts']:
                d=Path(attempt['directory'])
                if (d/'output_sha256.json').exists():rows.append(analyze_run(d,c,frozen['physics']))
    main=[r for r in rows if r['backend']=='HemoCell' and r['role']=='main' and r.get('accuracy') and r['accuracy'].get('measured_profile_si')]
    coarse=next((r for r in main if r['candidate_id']=='lbm_N16' and r['rank_count']==1),None)
    fine=next((r for r in main if r['candidate_id']=='lbm_N32'),None)
    refinement=None
    if coarse and fine:
        refinement=float(np.linalg.norm(np.array(coarse['accuracy']['measured_profile_si'])-fine['accuracy']['measured_profile_si'])/np.linalg.norm(fine['accuracy']['target_profile_si']))
    for r in rows:
        if r['role']!='main' or not r.get('accuracy'):continue
        a=r['accuracy'];ok=a.get('local_mean_flow_pass',False)
        if r['backend']=='HemoCell':
            passed=refinement is not None and refinement<=c['criteria']['refinement_relative_l2'] and coarse['accuracy']['local_mean_flow_pass'] and fine['accuracy']['local_mean_flow_pass']
            a['refinement_relative_l2']=refinement;a['gates']['refinement']=bool(passed)
            a['discretization_error_status']='TWO_GRID_CHECK_PASSED' if passed else 'INCONCLUSIVE_OR_FAILED'
            r['qualified']=bool(ok and passed)
            validation=fine if r['candidate_id']=='lbm_N16' else coarse
            r['qualification_validation_wall_s']=validation['timing']['total_wall_s'] if validation else None
        else:r['qualified']=bool(ok)
        if r['qualified']:
            r['time_to_qualified_solution_s']=r['timing']['total_wall_s']+(r.get('qualification_validation_wall_s') or 0)
            r['timing']['time_to_qualified_solution_s']=r['time_to_qualified_solution_s']
            r['qualification_time_scope']='Completed fixed window from fresh state plus other-grid validation once for LBM. Not an optimized earliest stopping time.'
    if not any(x['backend']=='Mirheo' for x in rows):
        g=frozen['gpu_plan'];rows.append(dict(backend='Mirheo',method='SDPD',case_id=c['case_id'],candidate_id='sdpd_target_linear',task_id='sdpd_main',role='main',
            source_category='NOT_RUN',source_commit=g['source_commit'],binary_hash=g['binary_hash'],precision=g['precision'],hardware='RTX 4060',rank_count=2,
            physical_parameters=frozen['physics'],dt_si=g['dt_si'],actual_steps=None,actual_time_si=None,timing=None,memory=None,accuracy=None,
            sampling_status='NOT_RUN',budget_status='PENDING_SCOPE_APPROVAL',qualified=False,time_to_qualified_solution_s=None,limitations=['Missing newly authorized matched trajectory; old data cannot supply a qualified speed ratio.']))
    speedups=[]
    mir=next(r for r in reversed(rows) if r['backend']=='Mirheo')
    for r in main:speedups.append({'task_id':r['task_id'],'qualified_speedup':qualified_ratio(mir,r)})
    return {'schema_version':1,'status':'COMPLETE_LOCAL_SCREEN' if mir['qualified'] and any(r['qualified'] for r in main) else 'PARTIAL',
            'case_id':c['case_id'],'results':rows,'refinement_relative_l2':refinement,'qualified_speedups':speedups,
            'human_review':'PENDING','pipe_benchmark':'PIPE_BENCHMARK_DEFERRED','pipe_reason':'Current SDPD wall kernel truncation, frozen-wall density and near-wall slip have not been validated; native DPD wall examples do not establish SDPD validity.',
            'scope':'Actual WSL deployment, CPU double LBM vs GPU single SDPD; not equal-device/precision algorithm ranking or RBC performance.'}


def historical_reference():
    path=PROJECT_ROOT/'data/sdpd_diagnostics/result_4ad869c67164f241/corrected_tasks.json'
    r=next(x for x in read_json(path) if x['task_id']=='sdpd_flow');p=r['performance']
    return {'source_category':'HISTORICAL_REFERENCE','source_path':str(path),'source_sha256':sha256_file(path),'task_id':'sdpd_flow',
            'actual_steps':p['actual_steps'],'actual_time_si':p['actual_physical_duration_s'],'total_wall_s':p['total_wall_s'],
            'compute_ms_per_step':p['synchronized_steady_compute_s_per_step']*1000,'workflow_wall_s_per_us':p['wall_s_per_physical_s']*1e-6,
            'temperature_mean_K':r['temperature']['mean']*298.15,'temperature_CI':None,'nu_relative_error':r['viscosity']['target_relative_error'],
            'profile_relative_l2':float(np.linalg.norm(np.array(r['viscosity']['measured_profile_star'])-r['viscosity']['target_profile_star'])/np.linalg.norm(r['viscosity']['target_profile_star'])),
            'qualified':False,'qualified_speedup':None,'limitations':'Different forcing, initialization stages, bins, outputs and sampling windows. Historical cost context only; never pooled into the matched case or a speedup.'}

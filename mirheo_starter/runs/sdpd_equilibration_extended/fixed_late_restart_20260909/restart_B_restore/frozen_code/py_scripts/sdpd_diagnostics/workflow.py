"""Source-bound CPU diagnosis and immutable report packages."""
import csv
import json
import math
import time
from pathlib import Path
import numpy as np
import yaml
from py_scripts.fluid_physics.common import (PROJECT_ROOT,read_json,write_json,sha256_file,
    fingerprint,output,resource_snapshot,environment_identity,now)
from py_scripts.fluid_physics.units import Units
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.runner import shared_budget_state
from py_scripts.fluid_comparison.experiments import load_config as load_comparison,verify_historical
from py_scripts.fluid_comparison.models import frozen_targets
from py_scripts.fluid_comparison.analysis import analyze_run
from py_scripts.fluid_comparison.reporting import assess
from py_scripts.fluid_physics.analysis import read_csv
from .calculations import temperature_audit,profile_audit,sampling_row,matched_flow_audit
from .native import source_contract,snapshot_audit,pressure_audit


def load_config(path):
    path=Path(path);path=path if path.is_absolute() else PROJECT_ROOT/path
    c=yaml.safe_load(path.read_text());c['_config_path']=str(path);c['_config_sha256']=sha256_file(path)
    package=Path(c['comparison_result']);verify_package(package)
    if sha256_file(package/'package_sha256.json')!=c['comparison_manifest_sha256']:raise ValueError('COMPARISON_PACKAGE_IDENTITY_MISMATCH')
    reproduction=Path(c['setup_directory'])/'reproduction_before_repairs.json'
    if sha256_file(reproduction)!=c['reproduction_sha256']:raise ValueError('REPRODUCTION_IDENTITY_MISMATCH')
    if c['policy']['original_total_GPU_budget_s']!=3600 or not c['policy']['no_new_budget']:raise ValueError('BUDGET_AUTHORIZATION_CHANGED')
    for key in ('output_root','review_root','runs_root'):c[key]=str(output(c[key]))
    return c


def scalar_metrics(t):
    v=t.get('viscosity') or {}
    return {'temperature_mean':t.get('temperature',{}).get('mean'),'temperature_CI_halfwidth':t.get('temperature',{}).get('ci95_halfwidth'),
        'pressure_mean_star':t.get('pressure',{}).get('mean'),'pressure_CI_halfwidth':t.get('pressure',{}).get('ci95_halfwidth'),
        'nu_si':v.get('nu_si'),'nu_CI95_si':v.get('ci95_si'),'profile_relative_RMS':v.get('relative_rms'),
        'viscosity_blocks':v.get('block_valid_count'),'stationarity_status':t.get('stationarity_status')}


def plan_experiments(c,rows,budget,u,tol):
    by={r['task_id']:r for r in rows};half=by['sdpd_half_dt'];base=by['sdpd_equilibrium'];cfg=c['probe']
    dt=base['parameters']['dt_star']*cfg['dt_factor'];steps=round(cfg['duration_star']/dt)
    speed=half['performance']['total_wall_s_per_step'];predicted=steps*speed*cfg['cost_safety_factor']+5
    usable=min(cfg['allocation_s'],budget['remaining_s'])-12
    affordable=predicted<=usable
    one={'id':cfg['id'],'status':'ELIGIBLE_WITHIN_EXISTING_AUTHORIZATION' if affordable else 'REQUIRES_ADDITIONAL_AUTHORIZATION',
       'unique_main_change':'Unforced dt is halved; same native model, units, Uniform constructor, mass, N, velocity seed and physical sample interval.',
       'baseline_task':'sdpd_equilibrium','baseline_reused_steps':round(cfg['duration_star']/base['parameters']['dt_star']),
       'dt_star':dt,'steps':steps,'physical_duration_s':cfg['duration_star']*u.t0,'predicted_total_wall_s':predicted,
       'measured_cost_source':'sdpd_half_dt total charged wall / actual steps','allocation_s':cfg['allocation_s'],
       'purpose':cfg['purpose'],'interpretation_rule':cfg['interpretation_rule'],'initial_state_rule':cfg['initial_state_rule'],
       'comparison_windows_star':cfg['predeclared_comparison_windows_star'],'stop_conditions':cfg['stop_conditions'],
       'qualification_capable':False,'why_still_informative':'It spans the observed heating peak and substantial decay, so it tests transient time-step dependence without repeating an inadequate qualification attempt.'}
    # Existing long-record ACF contains drift: these estimates are conditional,
    # not a claim that the needed stationary correlation time has been measured.
    flow=by['sdpd_flow'];stats=flow['viscosity']['slope_statistics'];block=stats['block_duration_star']
    floor=tol['block_relaxation_multiples']*flow['relaxation_prior_star']
    warmup=.30;estimate=[]
    flow_dt=flow['parameters']['dt_star'];every=flow['parameters']['snapshot_every'];spacing=flow_dt*every
    warmup_steps=int(math.ceil(warmup/spacing-1e-9))*every
    for label,b in [('optimistic_relaxation_floor',floor),('observed_drifting_ACF_upper_scenario',block)]:
        block_samples=int(math.ceil(b/spacing-1e-9));step_count=warmup_steps+8*block_samples*every
        wall=step_count*flow['performance']['total_wall_s_per_step']*1.15+18
        estimate.append({'scenario':label,'warmup_star':warmup_steps*flow_dt,'block_star':block_samples*spacing,
          'unrounded_block_requirement_star':b,'block_samples':block_samples,'sample_interval_star':spacing,'blocks':8,'steps':step_count,
          'duration_star':step_count*flow['parameters']['dt_star'],'physical_duration_s':step_count*flow['parameters']['dt_star']*u.t0,
          'predicted_wall_s':wall,'single_task_600s_feasible':wall<=600,
          'status':'REQUIRES_ADDITIONAL_AUTHORIZATION','additional_budget_needed_s':max(0.,wall-budget['remaining_s']),
          'assumptions':'Warmup of 0.30 star is a prospective design assumption; stationary correlation time and actual plateau have not been established. The ACF scenario includes trend and can overestimate stationary correlation.'})
    next_steps=400000;next_wall=next_steps*base['performance']['total_wall_s_per_step']*1.15+18
    priority={'id':'longer_unforced_thermal_plateau','status':'REQUIRES_ADDITIONAL_AUTHORIZATION','steps':next_steps,
      'dt_star':base['parameters']['dt_star'],'duration_star':next_steps*base['parameters']['dt_star'],
      'physical_duration_s':next_steps*base['parameters']['dt_star']*u.t0,'predicted_total_wall_s':next_wall,
      'recommended_allocation_s':600,'additional_authorization_needed_for_allocation_s':max(0.,600-budget['remaining_s']),
      'single_task_limit_s':600,'unique_main_change':'Longer unforced observation with unchanged native model, initialization, units and dt.',
      'prospective_statistics_start_star':.30,'prospective_observation_duration_star':.10,
      'predeclared_decision':'Use the entire predeclared final 0.10 star interval, require at least 8 blocks with block >= max(5 ACF, 2 viscous relaxation times), two-half thermal drift <=0.005 and |mean/target-1| + CI_halfwidth/target <=0.02. If these fail, record the failure; do not shift the window or extend automatically.',
      'assumptions':'The proposed 0.30-star warmup and stationary ACF are unverified. This run can establish or refute a plateau in its recorded interval; 8 blocks and qualification are not guaranteed. The 600-s allocation includes the original 12-s shutdown reserve.',
      'cost_source':'Old unforced sdpd_equilibrium total charged wall/actual steps, 15% safety factor, plus 18 seconds.'}
    return {'selected_probe':one,'other_probes':[
      {'id':'official_low_parameter_smoke','status':'NOT_SELECTED','reason':'The actual target model already ran; a low-parameter smoke cannot explain this transient.'},
      {'id':'pre_equilibrated_restart','status':'NOT_TESTED','reason':'Old files do not contain full velocities and random state; a valid continuation has not been established.'},
      {'id':'continuous_vs_chunked_run','status':'NOT_SELECTED','reason':'Pinned call-chain review finds no particle/time/SDPD RNG reset in the present pure-liquid setup. Statistical floating-point trajectory equivalence is not claimed.'},
      {'id':'low_output_vs_current','status':'NOT_SELECTED','reason':'Getters/savers are read-only to physical channels; CPU moment identities and synchronization checks address the immediate observability question.'},
      {'id':'box_size_or_smooth_force','status':'REQUIRES_ADDITIONAL_AUTHORIZATION','reason':'Defer until temperature and stationarity are established; finite rc/channel remains an unconfirmed residual cause.'}],
      'eight_block_cost_scenarios':estimate,'next_priority_experiment':priority,'no_automatic_extra_budget':True,
      'next_priority':'Establish a stable unforced target-liquid thermal plateau with sufficient independent blocks; any run exceeding 600 s first requires a verified state/random-process checkpoint strategy.'}


def build_evidence(temperatures,profiles,sampling,pressure,plan):
    by={t['task_id']:t for t in temperatures};flow=by['sdpd_flow'];eq=by['sdpd_equilibrium']
    s=next(x for x in sampling if x['task_id']=='sdpd_flow' and x['observable']=='viscosity_slope')
    rows=[]
    def add(id,phenomenon,state,explanation,support,against,missing,test,confidence='high'):
        rows.append({'id':id,'phenomenon':phenomenon,'status':state,'one_sentence':explanation,'supporting_evidence':support,
          'counter_evidence_or_limit':against,'missing_evidence':missing,'minimal_test':test,'confidence':confidence})
    add('temperature_moment_formula','温度测量','RULED_OUT_WITHIN_TESTED_SCOPE','保存的分箱矩可独立重算 worker 的原始、COM 和粗分箱温度，未发现能解释大幅升温的质量或自由度错误。',
        {'task':'sdpd_equilibrium','CPU_recovery':eq['independent_moment_recovery'],'N_minus_one_relative_effect':eq['DOF_relative_correction']},
        '原生 Stats 与 worker 相差一个积分步；二者读取同一底层速度，CPU 矩复核也不是另一套测量仪器。','旧逐粒子速度与分量方差未保存。','新短探针保存实际前后相位速度，验证 getter/通道和 kick 更新。')
    add('thermal_transient','温度偏高','CONFIRMED','选入统计的轨迹仍包含初始化后的降温过程，所报均值不能当作已证明稳定的平台温度。',
        {'equilibrium_peak':eq['peak_kBT'],'equilibrium_last_10pct':eq['last_10pct_mean'],'all_windows':eq['windows'],'flow_windows':flow['windows']},
        '这证实存在过渡段，尚不独立证明其全部由位置初始化而非积分偏差造成。','稳态 dt/dt/2 温度与重复初始化。','先比较相同物理时间内无驱动 dt/2 的升温/衰减趋势。')
    add('initial_structure','初始化结构与背景压力','SUPPORTED_BUT_NOT_CONFIRMED','随机位置、较大正背景配对压力与持续变化的核密度支持结构松弛驱动启动升温，但完整能量收支未保存。',
        {'source':'native_source_contract.json','snapshots':'snapshot_audit.json','series':'temperature_audit.json'},
        '固定快照的保守力和势能是公式推导量，不是已保存的全部原生力或连续能量轨迹。','t=0 力、耗散/随机功及充分平衡初态。','本轮先以无驱动半步趋势及新保存的前后状态缩小范围。','medium')
    add('integrator_bias','时间推进','SUPPORTED_BUT_NOT_CONFIRMED','本地纯液体调用链实际采用一次 kick-drift；固定耗散模态会有正温度步长偏差，但不能把简化估计直接套到液体。',
        {'source':'vv.cu / integration_kernel.h / Simulation scheduler','analytic_reference':'T_discrete/T = 1/(1-lambda*dt/2)','scales':'snapshot_audit.json'},
        '没有证据把原生库定为错误；相同种子、不同 dt 也不表示相同连续噪声轨迹。','稳定状态下的物理时长匹配 dt 对照。','预登记的无驱动半步趋势探针仅判断启动过渡的步长依赖。','medium')
    add('flow_subtraction','流动温度扣除','RULED_OUT_WITHIN_TESTED_SCOPE','平衡任务也升温，粗细分箱差异很小；漏扣平均 Poiseuille 流速不足以解释当前整体温升。',
        {'equilibrium_has_no_force':True,'flow_bin_empty_count':flow['empty_bin_records'],'flow_moments':flow['independent_moment_recovery']},
        '仍不能排除分箱内集体运动与未保存的局部速度异常。','逐粒子旧速度和独立空间/分量温度。','复核空间粗分箱矩；新快照只能补新运行的证据。')
    add('one_block','只有一个统计块','CONFIRMED','预算缩短轨迹、初始化截段及漂移信号的长 ACF 联合留下一个完整块；块数整除本身没有算错。',
        s,'漂移下估计的 ACF 不是可靠稳态相关时间，不能机械将当前总时长乘八。','真正稳定后的相关时间。','sampling_audit.csv 给出每项取整和尾部账单；proposed_experiments.json 给条件预算。')
    add('tail_estimator','均值与置信区间','CONFIRMED','旧均值包含不完整尾部，而区间只用完整块；已通过先失败后通过的回归测试修复，并保留全窗口均值。',
        {'test':'StatisticsRepairTests.test_mean_and_ci_use_same_complete_blocks','diff':'repair_differences.json'},
        '这只修正估计量的一致性，不消除过渡、相关性或物性误差。','无。','用完全相同原始数据重新计算；无需 GPU。')
    add('warmup_gate','平衡段判据','CONFIRMED','现有截段只搜索 20–50% 且允许 10% 两半漂移，这不足以证明满足 2% 温度要求的稳态。',
        {'windows':flow['windows']},'材料温度总门槛仍会拒绝候选，未将其最终标为合格。','事前定义的更严格稳态判据和更长观测。','报告所有窗口；本轮不修改阈值使历史记录通过。')
    add('profile_residual','黏度与剖面不符','SUPPORTED_BUT_NOT_CONFIRMED','拟合系数只描述一个投影；过渡和不足统计使剖面剩余形状无法可靠归因，有限作用范围也尚未排除。',
        {'profiles':'profile_audit.json','rc_over_half_channel':profiles[0]['rc_over_half_channel'] if profiles else None},
        '不能把残差位置特征直接证明为核截断或有限盒效应；不能事后删除变号附近点。','稳定多块的空间误差条，以及尺寸或平滑驱动对照。','先取得可靠热平衡，再做有限尺度诊断。','medium')
    add('pressure_coverage','压力覆盖','CONFIRMED','端点减基准的名义余量小于其统计不确定性；参考压力共享产生相关误差，移动零点不会增加覆盖能力。',
        {'margins':pressure['endpoint_margins'],'rounding_pa':pressure['CSV_virial_rounding_max_half_quantum_pa']},
        '各 EOS 点温度/结构还在松弛，当前区间也不等于合格等温 EOS；独立运行协方差假设只是诊断。','合格同温度下的 EOS 重复测量及运行间协方差。','先热平衡；随后按端点与基准的压差裕量设计精度。')
    add('repeated_run','重复 u.run','RULED_OUT_WITHIN_TESTED_SCOPE','本地调用链保留物理状态、累计时间和 SDPD 随机发生器；未发现每次采样重新初始化液体。',
        {'call_chain':'native_source_contract.json'},'最终 cell-list 重排可能改变浮点轨迹，未声称逐位等同连续 run。','相同初态连续/分块的实测对照未做。','目前源码证据足以排除状态重置，不为此消耗有限 GPU 预算。')
    add('missing_old_state','可观测性限制','INSUFFICIENT_EVIDENCE','旧位置快照不足以恢复逐粒子速度、力、随机状态或可靠续跑。',
        {'files':'provenance.json run manifests'},'矩恒等式可核查聚合温度，但不能还原粒子尾部分布或所有速度分量。','旧未保存的通道。','仅在新探针保存需要的真实状态；不伪造旧数据。')
    return rows


def export(c):
    cc=load_comparison(c['comparison_config']);case,mapping,u=frozen_targets(cc);verify_historical(cc)
    original=read_json(Path(c['comparison_result'])/'comparison_summary.json')
    reproduction=read_json(Path(c['setup_directory'])/'reproduction_before_repairs.json')
    if reproduction['status']!='PASS':raise ValueError('ORIGINAL_REPORT_NOT_REPRODUCED')
    budget=shared_budget_state(Path(cc['runs_root'])/cc['campaign_id'],cc);budget.pop('pool',None)
    native=source_contract();native_hashes={x['path']:x['sha256'] for x in native['sources']}
    roots=[PROJECT_ROOT/'py_scripts/sdpd_diagnostics',PROJECT_ROOT/'py_scripts/fluid_physics',PROJECT_ROOT/'py_scripts/fluid_comparison']
    code={str(p):sha256_file(p) for folder in roots for p in sorted(folder.glob('*.py'))}
    inputs={t['directory']:sha256_file(Path(t['directory'])/'output_sha256.json') for t in original['tasks']}
    for name in c['probe'].get('preserved_failed_probe_ids',[]):
        path=Path(c['runs_root'])/c['campaign_id']/name
        inputs[str(path)]=sha256_file(path/'output_sha256.json')
    probe_path=Path(c['runs_root'])/c['campaign_id']/c['probe']['id'];probe_manifest=probe_path/'output_sha256.json'
    identity={'config_sha256':c['_config_sha256'],'comparison_manifest':c['comparison_manifest_sha256'],'code':code,
              'native_sources':native_hashes,'raw_inputs':inputs,'budget':budget,
              'probe_manifest':sha256_file(probe_manifest) if probe_manifest.exists() else None}
    validation=Path(c['setup_directory'])/'CPU_validation.json'
    if validation.exists():identity['CPU_validation']=read_json(validation)
    d=output(Path(c['output_root'])/('result_'+fingerprint(identity)[:16]))
    if d.exists():verify_package(d);return d
    started=time.monotonic();oldc=yaml.safe_load((Path(cc['historical_result'])/'fluid_physics.yaml').read_text())
    oldc['legacy_windows']=cc['legacy_windows'];oldc['proposed_additional']=cc['proposed_additional']
    corrected=[];temperatures=[];profiles=[];sampling=[];differences=[];provenance=[]
    for old in original['tasks']:
        hist=old['source_category']=='HISTORICAL_REFERENCE';config=oldc if hist else cc
        fresh=analyze_run(old['directory'],config,u,case,historical=hist);corrected.append(fresh)
        spec=fresh['parameters'];execution=fresh['execution']
        provenance.append({'task_id':old['task_id'],'directory':old['directory'],'worker_sha256':sha256_file(Path(old['directory'])/'gpu_worker.py'),
          'parameters_sha256':sha256_file(Path(old['directory'])/'actual_parameters.json'),'library_sha256':spec.get('libmirheo_sha256'),
          'desired_steps':spec.get('desired_steps',spec['steps']),'allocated_steps':spec['steps'],'actual_steps':fresh['performance']['actual_steps'],
          'dt_star':spec['dt_star'],'snapshot_every':spec['snapshot_every'],'N':fresh['actual_N'],'actual_duration_s':fresh['performance']['actual_physical_duration_s'],
          'status':execution['status'],'timeout':execution.get('timeout'),'exit_code':execution.get('exit_code'),
          'worker_completion':read_json(Path(old['directory'])/'worker_completion.json'),'raw_manifest_sha256':inputs[old['directory']]})
        differences.append({'task_id':old['task_id'],'before':scalar_metrics(old),'after':scalar_metrics(fresh),'reason':'Analysis-only complete-block estimator alignment. Original trajectories and acceptance records remain immutable.'})
        if 'relaxation_prior_star' in fresh:
            for name,stats in [('temperature',fresh['temperature']),('pressure',fresh['pressure']),('viscosity_slope',fresh.get('viscosity',{}).get('slope_statistics'))]:
                if stats:
                    row=sampling_row(old['task_id'],spec,name,stats,fresh['relaxation_prior_star'],config['proposed_acceptance'],u,execution['elapsed_monotonic_s'])
                    row.update(actual_steps=fresh['performance']['actual_steps'],actual_physical_duration_s=fresh['performance']['actual_physical_duration_s'],
                               actual_duration_star=fresh['performance']['actual_steps']*spec['dt_star'],execution_status=execution['status'])
                    sampling.append(row)
        if old['method']=='SDPD' or old['task_id'] in ('equilibrium','high_gamma_equilibrium'):
            temperatures.append(temperature_audit(old['directory'],spec,old,fresh,u,c['window_fractions'],config['proposed_acceptance']))
        if old['method']=='SDPD' and old['kind']=='flow':profiles.append(profile_audit(old['directory'],spec,old,fresh,u))
    assessments=assess(corrected,case,u,cc);pressure=pressure_audit(original['candidate_assessments'],assessments,corrected,case,u)
    snapshots=[]
    for task_id in c['snapshot_tasks']:
        t=next(r for r in corrected if r['task_id']==task_id);paths=sorted(Path(t['directory']).glob('snapshot_*_positions.npy'))
        for p in [paths[0],paths[-1]]:snapshots.append(snapshot_audit(p,t['parameters']))
    plan=plan_experiments(c,corrected,budget,u,cc['proposed_acceptance'])
    evidence=build_evidence(temperatures,profiles,sampling,pressure,plan)
    matches=[]
    base=next(t for t in corrected if t['task_id']=='sdpd_flow')
    for name in ('sdpd_half_dt','sdpd_half_force'):
        other=next(t for t in corrected if t['task_id']==name)
        matches.append(matched_flow_audit(base,other,u))
    probe=None
    if probe_manifest.exists():
        from .probes import analyze_probe
        probe=analyze_probe(probe_path,c,corrected,u)
        plan['selected_probe']['status']=probe['summary']['status']
        plan['selected_probe']['actual_charged_s']=probe['summary']['charged_s']
        for e in evidence:
            if e['id'] in ('initial_structure','integrator_bias','temperature_moment_formula'):
                e['supporting_evidence']['new_probe']='diagnosis_summary.json / new_probe'
            if e['id']=='initial_structure' and probe['summary']['status']=='COMPLETED_FIXED_PHYSICAL_WINDOW':
                e['one_sentence']='真实初态温度约为目标值，dt 与 dt/2 均出现近似的强温升与衰减；结合核密度和保守势能下降，更支持初始结构松弛，完整因果能量收支尚未测得。'
                e['supporting_evidence']['peak_comparison']=[probe['summary']['baseline_peak'],probe['summary']['half_dt_peak']]
            if e['id']=='integrator_bias' and probe['summary']['status']=='COMPLETED_FIXED_PHYSICAL_WINDOW':
                e['counter_evidence_or_limit']='固定物理时段减半 dt 后强温升仍在，反对将启动峰值主要归于大步长；短探针没有检验稳态温偏。原生库缺陷仍未证实。'
    failed=[]
    for name in c['probe'].get('preserved_failed_probe_ids',[]):
        path=Path(c['runs_root'])/c['campaign_id']/name
        for relative,h in read_json(path/'output_sha256.json').items():
            if sha256_file(path/relative)!=h:raise ValueError('FAILED_ATTEMPT_HASH_MISMATCH')
        failed.append({'directory':str(path),'execution':read_json(path/'execution.json'),
          'raw_manifest_sha256':sha256_file(path/'output_sha256.json'),
          'cause':'New diagnostic saver requested forces; pinned native name is __forces. Failed before integration; no physical evidence.',
          'correction':'Diagnostic observer only; failing source-contract test preceded the minimal one-token fix.'})
    env=environment_identity()
    oldenv=read_json(Path(c['comparison_result'])/'provenance.json')['environment']
    for key in ('library_sha256','mirheo_commit','precision'):
        if env[key]!=oldenv[key]:raise ValueError('NATIVE_ENVIRONMENT_CHANGED '+key)
    timescales=[]
    for r in corrected:
        if r['method']!='SDPD':continue
        spec=r['parameters'];spacing=r['actual_n_star']**(-1/3);dt=spec['dt_star'];cs=spec['candidate']['sound_speed_star'];nu=spec['candidate']['nu_prior_star']
        m=read_csv(Path(r['directory'])/'moments.csv');max_speed=float(np.max(m['max_speed_star']))
        timescales.append({'task_id':r['task_id'],'particle_spacing_star':spacing,'acoustic_spacing_time_star':spacing/cs,
          'acoustic_box_crossing_time_star':spec['domain_star'][1]/cs,'viscous_spacing_time_star':spacing**2/nu,
          'box_viscous_decay_time_star':spec['domain_star'][1]**2/(4*math.pi**2*nu),
          'dt_over_acoustic_spacing_time':dt*cs/spacing,'dt_over_viscous_spacing_time':dt*nu/spacing**2,
          'max_observed_speed_star':max_speed,'max_single_step_displacement_over_spacing':dt*max_speed/spacing,
          'structure_relaxation_time_star':None,'structure_reason':'Snapshots and temperature still evolve; no verified stationary structural relaxation time.'})
    summary={'status':'DIAGNOSED_WITH_LIMITATIONS','selection':None,'human_review':'PENDING','original_reproduction':'PASS',
      'analysis_repair':'CONFIRMED_AND_CPU_TESTED','no_numerical_model_repair_claim':True,'shared_budget':budget,
      'corrected_candidate_assessments':assessments,'evidence':evidence,'new_probe':probe,
      'preserved_failed_diagnostic_attempts':failed,
      'CPU_wall_s':time.monotonic()-started,'units':mapping,'physical_case':case,
      'execution_policy':'No GPU starts from analyze-only, no default solver/parameter selection, no new budget.'}
    artifacts={'diagnosis_summary':summary,'evidence_matrix':evidence,'temperature_audit':temperatures,'pressure_audit':pressure,'time_scale_audit':timescales,
      'profile_audit':profiles,'snapshot_audit':snapshots,'proposed_experiments':plan,'repair_differences':differences,
      'matched_physical_time_audit':matches,'native_source_contract':native,'corrected_tasks':corrected,
      'provenance':{'created_at':now(),'identity':identity,'tasks':provenance,'environment':env,'original_result':c['comparison_result'],'setup':c['setup_directory']},
      'resource_record':resource_snapshot(),'unit_audit':{'density_si':u.to_si(mapping['rho_star'],'mass_density'),
       'temperature_K':u.temperature_K(mapping['kBT_star']),'target_nu_star':mapping['required_nu_star'],'target_mu_star':mapping['required_mu_star'],
       'mu_over_nu_star':mapping['required_mu_star']/mapping['required_nu_star'],'mass_density_star':mapping['rho_star'],
       'scales_si_per_star':u.scales,'particle_mass_si':mapping['m_star']*u.M0,'system_mass_si_at_N4096':4096*mapping['m_star']*u.M0,
       'baseline_flow_force_si':float(u.to_si(base['parameters']['task']['force_star'],'force')),
       'baseline_flow_acceleration_si':float(u.to_si(base['parameters']['task']['force_star']/base['parameters']['m_star'],'acceleration')),
       'temperature_algebra_is_not_qualification':True}}
    # Validate complete payloads before creating the immutable output directory.
    # Convert NumPy scalar representations only; never replace invalid physics by zeros.
    def plain(value):
        if isinstance(value,np.generic):return value.item()
        if isinstance(value,np.ndarray):return plain(value.tolist())
        if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [plain(v) for v in value]
        return value
    artifacts=plain(artifacts)
    json.dumps(artifacts,allow_nan=False)
    d.mkdir(parents=True,exist_ok=False)
    for name,value in artifacts.items():write_json(d/(name+'.json'),value)
    for name,rows in [('sampling_audit',sampling),('temperature_audit',[{k:v for k,v in row.items() if isinstance(v,(str,int,float)) or v is None} for row in temperatures])]:
        columns=list(dict.fromkeys(k for row in rows for k in row))
        with (d/(name+'.csv')).open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    (d/'sdpd_diagnostics.yaml').write_text(Path(c['_config_path']).read_text())
    from .report import chinese_report,repair_log
    (d/'report_zh.md').write_text(chinese_report(summary,artifacts,sampling))
    (d/'repair_log.md').write_text(repair_log(differences,c))
    write_json(d/'package_sha256.json',{p.name:sha256_file(p) for p in d.iterdir() if p.is_file()});verify_package(d)
    return d

"""Prepare, authorize, run and reanalyze one fixed unforced SDPD experiment."""
import copy
import csv
import json
import math
from pathlib import Path
import shlex
import shutil
import numpy as np
import yaml
from py_scripts.fluid_physics.common import (PROJECT_ROOT, read_json, write_json, sha256_file,
    fingerprint, output, atomic_state, environment_identity, now)
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.runner import exclusive_lock, shared_budget_state, run_attempt
from py_scripts.fluid_comparison.experiments import load_config as comparison_config
from .equilibration_analysis import empty_result, analyze_series
from .equilibration_worker import RAW_COLUMNS, PRESSURE_DEFINITION, mechanical_pressure
from .observations import phase_audit


def load_config(path):
    p = Path(path).resolve()
    c = yaml.safe_load(p.read_text())
    if c['schema_version'] != 1 or c['extra_authorized_gpu_seconds'] != 0:
        raise ValueError('CONFIG_CANNOT_GRANT_GPU_AUTHORIZATION')
    for key in ['campaign_id','case_id']:
        if not c[key] or '/' in c[key] or c[key] in ['.','..']:
            raise ValueError('INVALID_RUN_ID')
    for key in ['diagnostic_config','delivery_record','diagnostic_package','data_root','runs_root','review_root','setup_directory','cpu_validation_record']:
        c[key] = str(output(c[key]))
    if c['policy']['selection'] is not None or c['policy']['single_task_limit_s'] != 600 or c['policy']['GPU_concurrency'] != 1:
        raise ValueError('INVALID_SCOPE_OR_BUDGET_POLICY')
    criteria=c['criteria']
    if (criteria['min_blocks']<8 or criteria['min_half_blocks']<4 or criteria['ACF_multiplier']!=5
            or criteria['block_relaxation_multiples']<2 or criteria['min_observation_duration_star']<.1
            or criteria['temperature_two_half_drift_relative']>.005
            or criteria['temperature_relative_error_including_CI']>.02):
        raise ValueError('CRITERIA_CANNOT_WEAKEN_PREREGISTERED_SCREEN')
    if any(not math.isfinite(v) or v<=0 for v in criteria.values() if isinstance(v,(int,float))):
        raise ValueError('INVALID_SCREENING_THRESHOLD')
    c['_config_path'], c['_config_sha256'] = str(p), sha256_file(p)
    return c


def code_identity():
    paths = sorted((PROJECT_ROOT/'py_scripts').rglob('*.py'))
    paths += [PROJECT_ROOT/'py_scripts/sdpd_equilibration.yaml']
    return {str(p):sha256_file(p) for p in paths if '__pycache__' not in p.parts}


def verify_manifest(directory, name='output_sha256.json'):
    d = Path(directory).resolve()
    hashes = read_json(d/name)
    for relative, digest in hashes.items():
        p = (d/relative).resolve()
        if not p.is_relative_to(d) or sha256_file(p) != digest:
            raise ValueError('IMMUTABLE_RAW_HASH_MISMATCH '+str(p))


def source_bundle(c):
    record_path, package = Path(c['delivery_record']), Path(c['diagnostic_package'])
    if sha256_file(record_path) != c['delivery_record_sha256']:
        raise ValueError('DIAGNOSTIC_DELIVERY_RECORD_CHANGED')
    delivery = read_json(record_path)
    if (delivery['status'] != 'PASS' or Path(delivery['package']).resolve() != package
            or delivery['package_manifest_sha256'] != c['diagnostic_manifest_sha256']
            or sha256_file(package/'package_sha256.json') != c['diagnostic_manifest_sha256']):
        raise ValueError('DELIVERED_PACKAGE_IDENTITY_MISMATCH')
    verify_package(package)
    diag = yaml.safe_load(Path(c['diagnostic_config']).read_text())
    provenance = read_json(package/'provenance.json')
    if sha256_file(Path(c['diagnostic_config'])) != provenance['identity']['config_sha256']:
        raise ValueError('HISTORICAL_DIAGNOSTIC_CONFIG_CHANGED')
    comparison = Path(diag['comparison_result'])
    verify_package(comparison)
    old = next(r for r in read_json(package/'corrected_tasks.json') if r['task_id'] == 'sdpd_equilibrium')
    raw = Path(old['directory'])
    verify_manifest(raw)
    spec = read_json(raw/'actual_parameters.json')
    if spec != old['parameters']:
        raise ValueError('HISTORICAL_ACTUAL_PARAMETERS_MISMATCH')
    proposal = read_json(package/'proposed_experiments.json')['next_priority_experiment']
    case, units = read_json(comparison/'physical_case.json'), read_json(comparison/'unit_mapping.json')
    env = environment_identity()
    for key in ['library_sha256','mirheo_commit','precision']:
        if env[key] != provenance['environment'][key]:
            raise ValueError('NATIVE_ENVIRONMENT_CHANGED '+key)
    for p, digest in provenance['environment']['sha256'].items():
        if sha256_file(Path(p)) != digest:
            raise ValueError('COMPATIBILITY_OR_NATIVE_INPUT_CHANGED '+p)
    for row in read_json(package/'native_source_contract.json')['sources']:
        if sha256_file(Path(row['path'])) != row['sha256']:
            raise ValueError('PINNED_NATIVE_CALL_CHAIN_CHANGED')
    return dict(spec=spec, proposal=proposal, case=case, units=units, environment=env,
                raw=raw, completion=read_json(raw/'worker_completion.json'), execution=read_json(raw/'execution.json'),
                comparison_config=comparison_config(diag['comparison_config']), provenance=provenance)


def exact_steps(duration, dt):
    steps = round(duration/dt)
    if steps <= 0 or not math.isclose(steps*dt,duration,abs_tol=1e-12,rel_tol=0):
        raise ValueError('DURATION_NOT_ON_STEP_GRID')
    return steps


def build_plan(c, source):
    p, old, proposal = c['plan'], source['spec'], source['proposal']
    if not (p['steps']==400000 and p['dt_star']==old['dt_star']==1e-6
            and p['formal_window_star']==[.3,.4] and p['duration_star']==.4 and p['allocation_s']==600
            and p['steps']==proposal['steps']==exact_steps(p['duration_star'],p['dt_star'])
            and p['dt_star']==proposal['dt_star']
            and math.isclose(p['sample_interval_star'],old['snapshot_every']*old['dt_star'],abs_tol=1e-15,rel_tol=0)):
        raise ValueError('MAIN_EXPERIMENT_DIFFERS_FROM_FROZEN_PROPOSAL')
    case, units = source['case'], source['units']
    if old['candidate']['viscosity_mu_star'] != units['required_mu_star']:
        raise ValueError('FROZEN_VISCOSITY_MAPPING_MISMATCH')
    for key in ['L0','M0','t0','E0']:
        if old['locked_units'][key] != units[key]:raise ValueError('FROZEN_UNIT_MISMATCH')
    spec = copy.deepcopy(old)
    spec['historical_template'] = {'task':spec.pop('task'),'original_allocated_steps':old['steps'],
                                 'original_desired_steps':old['desired_steps'],'original_directory':str(source['raw'])}
    spec['task'] = copy.deepcopy(old['task'])
    spec['task'].update(id=c['case_id'],kind='equilibrium',test_kind='unforced_equilibration',
                        force_star=0.,desired_time_star=p['duration_star'],desired_steps=p['steps'],allocation_s=p['allocation_s'])
    spec.update(steps=p['steps'], desired_steps=p['steps'], desired_time_star=p['duration_star'],
                native_integrator='VelocityVerlet', source_category='FIXED_UNFORCED_EQUILIBRATION',
                expected_N=source['completion']['actual_N'], target_temperature_K=case['temperature_K'],
                script_sha256=sha256_file(Path(__file__).with_name('equilibration_worker.py')),
                snapshot_marks_steps=[exact_steps(x,p['dt_star']) for x in p['snapshot_times_star']],
                formal_window_star=p['formal_window_star'])
    if spec['task']['force_star'] != 0 or old['task']['force_star'] != 0 or spec['candidate']['method'] != 'SDPD':
        raise ValueError('UNFORCED_NATIVE_SDPD_REQUIRED')
    completion, execution = source['completion'], source['execution']
    base_rate = execution['elapsed_monotonic_s']/completion['steps']
    raw_estimate = p['steps']*base_rate
    predicted = raw_estimate*p['cost_safety_factor']+p['extra_setup_output_exit_reserve_s']
    if not math.isclose(predicted,proposal['predicted_total_wall_s'],abs_tol=1e-8):
        raise ValueError('COST_SOURCE_DOES_NOT_REPRODUCE_DELIVERED_PROPOSAL')
    minimum = c['criteria']['block_relaxation_multiples']*spec['domain_star'][1]**2/(4*math.pi**2*spec['candidate']['nu_prior_star'])
    plan = {**p,'parameters':spec,'criteria':copy.deepcopy(c['criteria']),
            'minimum_block_duration_star':minimum,'minimum_block_samples':math.ceil(minimum/p['sample_interval_star']),
            'physical_duration_s':p['duration_star']*spec['locked_units']['t0'],
            'formal_window_physical_s':[x*spec['locked_units']['t0'] for x in p['formal_window_star']],
            'planned_sample_count_including_t0':p['steps']//spec['snapshot_every']+1,
            'planned_formal_samples':round((p['formal_window_star'][1]-p['formal_window_star'][0])/p['sample_interval_star']),
            'predicted_total_wall_s':predicted,'predicted_unmargined_wall_s':raw_estimate,
            'cost_evidence':{'directory':str(source['raw']),'actual_steps':completion['steps'],
                'charged_s':execution['elapsed_monotonic_s'],'seconds_per_step_including_overhead':base_rate,
                'setup_s':completion['setup_wall_s'],'compute_s':completion['compute_wall_s'],
                'sampling_s':completion['snapshot_wall_s'],
                'other_launch_exit_s':execution['elapsed_monotonic_s']-completion['setup_wall_s']-completion['compute_wall_s']-completion['snapshot_wall_s'],
                'safety_margin_s':raw_estimate*(p['cost_safety_factor']-1),
                'additional_reserve_s':p['extra_setup_output_exit_reserve_s'],
                'new_observer_cost_measured':False,
                'caveat':'Same cadence, one extra persistent device stress saver per step plus sampled stress/IDs; 15% + 18 s are planning margins, not a measured new runtime.'},
            'normal_exit_deadline_s':p['allocation_s']-p['graceful_margin_s'],
            'fits_task_deadline_estimate':predicted<=p['allocation_s']-p['graceful_margin_s'],
            'continuity':'One coordinator and one particle initialization; matching MPI run chunks, no state/RNG reset.',
            'restart':'NOT_VERIFIED_NOT_USED; old trajectories are not prepended.',
            'stop_rules':['NaN/Inf','illegal density','particle count/ID change','resource limits','budget watchdog','planned steps'],
            'temperature_excursion_is_stop_condition':False,'pressure_definition':PRESSURE_DEFINITION,
            'case_fluid_values':{k:case[k] for k in ['density_kg_m3','kinematic_viscosity_m2_s','dynamic_viscosity_pa_s','temperature_K']},
            'source_delivery_sha256':c['delivery_record_sha256'],'source_manifest_sha256':c['diagnostic_manifest_sha256'],
            'observation_native_source_sha256':{str(PROJECT_ROOT/'vendor/Mirheo/src/mirheo'/p):sha256_file(PROJECT_ROOT/'vendor/Mirheo/src/mirheo'/p) for p in
                ['core/interactions/pairwise/sdpd.cu','plugins/particle_channel_saver.cpp','core/datatypes.h','bindings/particle_vectors.cpp','core/simulation.cpp']},
            'code_sha256':code_identity(),'selection':None}
    plan['plan_sha256'] = fingerprint(plan)
    return plan


def budget_status(c, source):
    cc = copy.deepcopy(source['comparison_config'])
    cc['read_only_unregistered_scope'] = True
    state = shared_budget_state(Path(c['runs_root'])/c['campaign_id'],cc,task_id=c['case_id'])
    state['ledger_sha256'] = {str(Path(row['directory'])/'budget_ledger.json'):sha256_file(Path(row['directory'])/'budget_ledger.json')
                              for row in state['members'] if (Path(row['directory'])/'budget_ledger.json').exists()}
    state['pool_sha256'] = sha256_file(Path(state['pool_path']))
    return state


def request_budget(c, plan, state):
    needed = max(0.,plan['allocation_s']-state['remaining_s'])
    return {'status':'READY_TO_RUN_AWAITING_BUDGET' if needed else 'WITHIN_RECORDED_AUTHORIZATION',
            'expected_wall_s':plan['predicted_total_wall_s'],'suggested_task_limit_s':plan['allocation_s'],
            'original_authorization_s':state['original_limit_s'],'original_remaining_s':state['original_remaining_s'],
            'current_scope_remaining_s':state['remaining_s'],'total_charged_or_reserved_s':state['total_charged_or_reserved_s'],
            'extra_authorized_gpu_seconds':state['extra_authorized_gpu_seconds'],
            'additional_required_s':needed,'additional_request_whole_seconds':math.ceil(needed),
            'authorization_scope':{'campaign_directory':str(Path(c['runs_root'])/c['campaign_id']),
                                   'task_id':c['case_id'],'plan_sha256':plan['plan_sha256']},
            'single_task_limit_s':600,'GPU_concurrency':1,'automatic_extension':False,
            'actual_spending_may_be_lower':True,'user_approval_received':False if needed else None,
            'budget_snapshot':state}


def cpu_validation(c, required=False):
    p = Path(c['cpu_validation_record'])
    if not p.exists():
        if required:raise ValueError('CPU_VALIDATION_REQUIRED')
        return None
    record = read_json(p)
    if record['status'] != 'PASS':raise ValueError('CPU_VALIDATION_FAILED')
    for name, digest in record['code_sha256'].items():
        if sha256_file(Path(name)) != digest:raise ValueError('CPU_VALIDATION_STALE '+name)
    if record['code_sha256'] != code_identity():raise ValueError('CPU_VALIDATION_CODE_SET_CHANGED')
    return record


def read_raw(directory):
    with (Path(directory)/'raw_statistics.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:raise ValueError('RAW_STATISTICS_EMPTY')
    return {key:np.array([float(row[key]) if row[key] else np.nan for row in rows]) for key in RAW_COLUMNS}


def audit_snapshots(directory, plan):
    from .native import snapshot_audit
    d, s, c = Path(directory), plan['parameters'], plan['criteria']
    identities = read_json(d/'initial_state_identity.json')
    for name,h in identities['sha256'].items():
        if sha256_file(d/name)!=h:raise ValueError('INITIAL_STATE_HASH_MISMATCH')
    initial_ids = np.load(d/'initial_ids.npy',allow_pickle=False)
    frames = read_json(d/'equilibration_observer.json')['snapshots']
    if [x['step'] for x in frames] != s['snapshot_marks_steps']:
        raise ValueError('FIXED_SNAPSHOT_SEQUENCE_INCOMPLETE')
    audit = []
    raw=read_raw(d)
    for frame in frames:
        prefix = d/f"snapshot_{frame['index']:02d}"
        ids = np.load(str(prefix)+'_ids.npy',allow_pickle=False)
        if sha256_file(Path(str(prefix)+'_ids.npy')) != frame['ids_sha256'] or not np.array_equal(np.sort(ids),np.sort(initial_ids)):
            raise ValueError('SNAPSHOT_PARTICLE_ID_ERROR')
        a = [np.load(str(prefix)+'_'+name+'.npy',allow_pickle=False) for name in
             ['positions','pre_velocities','pre_forces','post_positions','post_velocities']]
        stress_path=Path(str(prefix)+'_stresses.npy')
        if sha256_file(stress_path)!=frame['stresses_sha256']:raise ValueError('SAVED_STRESS_HASH_MISMATCH')
        measured_pressure,_=mechanical_pressure(np.load(stress_path,allow_pickle=False),a[-1],s['m_star'],float(np.prod(s['domain_star'])),s['dt_star'],frame['step'])
        pressure_rows=raw['pressure_star'][raw['step']==frame['step']]
        if len(pressure_rows)!=1 or not math.isclose(measured_pressure,float(pressure_rows[0]),rel_tol=1e-12,abs_tol=1e-10):
            raise ValueError('SNAPSHOT_PRESSURE_REDUCTION_MISMATCH')
        phase = phase_audit(*a,s['dt_star'],s['m_star'],s['domain_star'])
        q = snapshot_audit(Path(str(prefix)+'_positions.npy'),s)
        if (q['density_relative_error_max'] > c['CPU_kernel_relative_tolerance']
                or phase['kick_max_abs'] > c['phase_kick_absolute_tolerance']
                or phase['drift_PBC_max_abs'] > c['phase_drift_absolute_tolerance']):
            raise ValueError('SAME_PHASE_NATIVE_STATE_CPU_AUDIT_FAILED')
        audit.append({'step':frame['step'],'phase':phase,'same_frame_pressure_star':measured_pressure,**q})
    late = [x for x in audit if x['step']*s['dt_star'] >= plan['formal_window_star'][0]-1e-12]
    if len(late) < 3:raise ValueError('FORMAL_STRUCTURE_SNAPSHOTS_MISSING')
    nearest = np.array([x['nearest_distance_quantiles_star'][1:4] for x in late])/s['rc_star']
    neighbors = np.array([x['neighbors_mean'] for x in late])/(s['task']['n_star']*4*math.pi*s['rc_star']**3/3)
    nearest_change = float(np.max(np.ptp(nearest,axis=0)))
    neighbor_change = float(np.ptp(neighbors))
    passed = (nearest_change <= c['snapshot_nearest_quantile_drift_over_rc']
              and neighbor_change <= c['snapshot_neighbor_drift_over_expected_neighbors'])
    return {'measurement_status':'PASS','status':'PASS' if passed else 'INCONCLUSIVE',
            'nearest_quantile_range_over_rc':nearest_change,'neighbor_mean_range_over_expected_neighbors':neighbor_change,
            'snapshots':audit,'scope':'Three fixed late structure snapshots plus the full sampled kernel distribution; not spatial independence or continuum convergence.'}


def report_text(plan, request, result, package):
    s, b = result['summary'], request
    measured = s['actual_steps'] is not None
    lines = ['# 无驱动 SDPD 稳定性实验', '',
        f"状态：{s['status']}。执行：{s['execution_status']}。",
        '程序完成、液体稳定、稳定温度正确、全部物性合格分别记录；selection=null，人工验收 PENDING。', '',
        f"固定 {plan['steps']} 步，dt*={plan['dt_star']}，连续 {plan['physical_duration_s']*1e6:.9f} µs；正式窗口 t*∈(0.30,0.40]。",
        '0.30* 是事前设计的观察起点，不是已证明的平衡时间。启动峰值完整保留，不按温度偏离 2% 停机。',
        f"原参数 μ*={plan['parameters']['candidate']['viscosity_mu_star']!r}，盒 [8,8,8]，N=4096，n*=8，m*=1，kBT*=1，rc*=1，Linear EOS，cs*=120，rho0*=0。",
        f"原单位 L0={plan['parameters']['locked_units']['L0']!r} m，M0={plan['parameters']['locked_units']['M0']!r} kg，t0={plan['parameters']['locked_units']['t0']!r} s；目标 298.15 K。",
        '初始化只做一次 Uniform、高斯速度及 COM 扣除；保留原生压力、随机和耗散作用，无驱动力、无追加恒温器、无持续速度缩放。', '',
        '## 预算与实测成本来源', '',
        f"旧无驱动任务：{plan['cost_evidence']['actual_steps']} 步，实收 {plan['cost_evidence']['charged_s']:.6f} s。",
        f"等比例成本 {plan['predicted_unmargined_wall_s']:.3f} s；加 15% 余量及 18 s 初始化/输出/退出余量，预计 {plan['predicted_total_wall_s']:.3f} s。",
        f"建议单任务上限 {plan['allocation_s']} s，包含 12 s 合作退出预留；预计正常运行必须在 {plan['normal_exit_deadline_s']} s 前完成。",
        f"原授权剩余 {b['original_remaining_s']:.9f} s；本任务可用 {b['current_scope_remaining_s']:.9f} s；需要追加 {b['additional_required_s']:.9f} s，建议申请整数 {b['additional_request_whole_seconds']} s。",
        '新观测增加少量应力读取、粒子 ID 与固定快照，尚无新实测速度；余量不是耗时或获得稳定平台的保证。预算不足则拒绝启动，不缩短主实验。', '',
        '## 事前判据', '',
        '先看整个固定窗口的两半与四个固定块趋势，再估计相关时间。明显持续漂移时不报告可信稳态 ACF。',
        f"块长至少 max(5×新窗口及两半的 ACF，{plan['minimum_block_duration_star']:.12g}*)；至少 8 个完整块，两半各至少 4 块。",
        '均值与 CI 使用同一批完整块，尾部单列。每块均值用点显示，整条轨迹均值的 CI 用带显示，不伪造单块独立 CI。',
        '温度两半均值漂移加误差上界不超过目标热能的 0.5%；固定四分窗范围不超过 1%。温度初筛：|均值/目标−1|+CI半宽/目标≤2%。',
        f"压力漂移尺度为 n*kBT={plan['parameters']['task']['n_star']*plan['parameters']['kBT_star']}*（{plan['parameters']['task']['n_star']*plan['parameters']['kBT_star']*plan['parameters']['locked_units']['si_per_star']['pressure']:.9f} Pa），不以背景压力作分母。",
        '核密度均值以原 n*、方差以 n*²、分位数以 n* 归一化；COM 漂移以 sqrt(kBT/m) 归一化。精确阈值见 run_plan.json。',
        '这些是本次事前筛选设计，不是普适物理定律。压力稳定不等于出口压力覆盖，N*m/V 恒定不证明局部结构稳定。', '',
        '## 本次证据', '', s['reason'],
        f"液体稳定：{s['liquid_stationary']}；稳定温度匹配：{s['stable_temperature_matches_target']}。",
        f"实际步数：{s['actual_steps']}；温度均值*：{s['temperature_mean_star']}；温度 CI95*：{s['temperature_CI95_star']}。",
        '未运行时所有新测量字段为 null；raw_statistics.csv 只有表头。历史曲线单独标注，不能填成新结果。' if not measured else
        '实际结果按原固定窗口分析；缺少样本或误差过大时保持不确定，不移动窗口、不自动延长。',
        '未新增黏度或 EOS 验证，不能据此宣称全部液体物性合格。', '',
        '## 唯一下一步', '',
        ('明确追加授权后执行这一项完整无驱动实验。' if not measured else
         ('在稳定状态下进行无驱动 dt/dt2 对照（另行授权）。' if s['status']=='STATIONARY_BUT_THERMALLY_BIASED' else
          '根据本次固定窗口证据安排下一项验证；不自动启动流动、密度扫描或参数搜索。')), '',
        '```bash','cd /home/lzy/projects/mirheo_starter',
        '.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration.yaml --execute','```', '',
        '## 可追溯与保护', '',
        f"诊断来源由交付记录及 manifest 哈希定位；本次计划哈希 `{plan['plan_sha256']}`。",
        '复用原 runner 的单 GPU 锁、monotonic 看门狗与本进程组清理。旧用量、新用量、未完成预约统一核算；原授权记录和历史账本不覆盖。',
        '只在用户明确给出秒数和范围后追加授权记录；配置 extra_authorized_gpu_seconds=0 不能自动放行。',
        '普通快照不含完整原生随机/积分状态，不称为 checkpoint，也不拼接旧短轨迹。',
        f"压力定义：{PRESSURE_DEFINITION}",
        f"数据包：{package}"]
    return '\n'.join(lines)+'\n'


def export(c, analyze=False):
    source = source_bundle(c)
    plan = build_plan(c,source)
    state = budget_status(c,source)
    request = request_budget(c,plan,state)
    validation = cpu_validation(c)
    directory = Path(c['runs_root'])/c['campaign_id']/c['case_id']
    raw_hash = None
    if directory.exists():
        if not (directory/'output_sha256.json').exists():raise ValueError('UNFINISHED_ATTEMPT_REQUIRES_MANUAL_INSPECTION_NO_RETRY')
        verify_manifest(directory)
        raw_hash = sha256_file(directory/'output_sha256.json')
        executed_plan = read_json(directory/'run_plan.json')
        # CPU reanalysis may fix analysis code without changing the registered physics/window.
        for key in ['parameters','criteria','steps','dt_star','formal_window_star','duration_star']:
            if executed_plan[key] != plan[key]:raise ValueError('EXECUTED_PLAN_CHANGED_NO_REINTERPRETATION')
        if read_json(directory/'actual_parameters.json') != executed_plan['parameters']:
            raise ValueError('ACTUAL_PARAMETERS_NOT_THE_REGISTERED_PLAN')
        execution = read_json(directory/'execution.json')
        completion = read_json(directory/'worker_completion.json') if (directory/'worker_completion.json').exists() else {}
        try:
            full = completion.get('steps') == plan['steps'] and completion.get('status') == 'COMPLETED_PLANNED_STEPS'
            structure = audit_snapshots(directory,plan) if full else {}
            data = read_raw(directory) if (directory/'raw_statistics.csv').exists() else {}
            result = analyze_series(data,plan,execution,completion,structure)
        except (ValueError,OSError,KeyError) as exc:
            result = empty_result('MEASUREMENT_OR_ANALYSIS_ERROR',str(exc))
            result['summary'].update(actual_steps=completion.get('steps'), reason=str(exc),
                analysis_recovery='Reanalyze immutable raw files after a reproducible regression; rerun only if required raw channels are absent.')
    else:
        result = empty_result(request['status'])
        if not validation:result['summary']['status']='PREPARATION_CPU_VALIDATION_REQUIRED'
    identity = {'plan_sha256':plan['plan_sha256'],'budget_sha256':fingerprint(request),
                'raw_manifest_sha256':raw_hash,'CPU_validation_sha256':fingerprint(validation) if validation else None}
    run_id = ('result_' if raw_hash else 'preflight_')+fingerprint(identity)[:16]
    package = output(Path(c['data_root'])/run_id)
    if package.exists():verify_package(package);return package
    package.mkdir(parents=True,exist_ok=False)
    provenance = {'created_at':now(),'identity':identity,'environment':source['environment'],
                  'delivery_record':c['delivery_record'],'delivery_sha256':c['delivery_record_sha256'],
                  'diagnostic_package':c['diagnostic_package'],'diagnostic_manifest_sha256':c['diagnostic_manifest_sha256'],
                  'source_actual_parameters':str(source['raw']/'actual_parameters.json'),
                  'source_actual_parameters_sha256':sha256_file(source['raw']/'actual_parameters.json'),
                  'new_raw_directory':str(directory) if raw_hash else None,'CPU_validation':validation,
                  'config_path':c['_config_path'],'config_sha256':c['_config_sha256'],
                  'selection':None,'human_review':'PENDING'}
    artifacts = {'run_plan':plan,'budget_request':request,'equilibration_summary':result['summary'],
                 'stationarity_checks':result['stationarity'],'sampling_audit':result['sampling'],
                 'provenance':provenance,'block_statistics':result['statistics']}
    for name,value in artifacts.items():write_json(package/(name+'.json'),value)
    if raw_hash:
        if (directory/'raw_statistics.csv').exists():shutil.copyfile(directory/'raw_statistics.csv',package/'raw_statistics.csv')
        else:
            with (package/'raw_statistics.csv').open('x') as f:csv.writer(f).writerow(RAW_COLUMNS)
        shutil.copyfile(directory/'actual_parameters.json',package/'actual_parameters.json')
        shutil.copyfile(directory/'execution.json',package/'execution.json')
    else:
        with (package/'raw_statistics.csv').open('x') as f:csv.writer(f).writerow(RAW_COLUMNS)
    (package/'report_zh.md').write_text(report_text(plan,request,result,package))
    shutil.copyfile(c['_config_path'],package/'sdpd_equilibration.yaml')
    write_json(package/'package_sha256.json',{p.name:sha256_file(p) for p in package.iterdir() if p.is_file()})
    verify_package(package)
    return package


def execute(c):
    package = export(c)
    plan, request = read_json(package/'run_plan.json'), read_json(package/'budget_request.json')
    campaign = Path(c['runs_root'])/c['campaign_id']
    if (campaign/c['case_id']).exists():
        print('EXISTING_ATTEMPT_REANALYZED_NO_GPU',flush=True)
        return package
    cpu_validation(c,required=True)
    if request['additional_required_s'] > 1e-9:
        raise PermissionError(f"READY_TO_RUN_AWAITING_BUDGET: remaining={request['current_scope_remaining_s']:.9f}s; required=600s; additional={request['additional_required_s']:.9f}s; NO_GPU_STARTED")
    if not plan['fits_task_deadline_estimate']:
        raise ValueError('PREDICTED_RUNTIME_EXCEEDS_TASK_LIMIT_NO_LAUNCH')
    source = source_bundle(c)
    cc = copy.deepcopy(source['comparison_config'])
    cc['campaign_id'] = c['campaign_id']
    pool_path = Path(cc['shared_budget_pool'])
    with exclusive_lock(read_json(pool_path)['gpu_lock']):
        state = budget_status(c,source)
        if state['remaining_s']+1e-9 < plan['allocation_s'] or any(m['running_reservation_s'] for m in state['members']):
            raise RuntimeError('SHARED_BALANCE_CHANGED_NO_LAUNCH')
        for record in state['authorization_records']:
            scope = record['scope']
            if scope['campaign_directory']==str(campaign) and scope['task_id']==c['case_id'] and scope['plan_sha256']!=plan['plan_sha256']:
                raise ValueError('APPROVAL_BINDS_A_DIFFERENT_PLAN')
        pool = state['pool']
        if not any(Path(m['directory']).resolve()==campaign.resolve() for m in pool['members']):
            campaign.mkdir(parents=True,exist_ok=True)
            write_json(campaign/'shared_pool_before_registration.json',pool)
            pool['members'].append({'directory':str(campaign),'campaign_id':c['campaign_id'],'ledger_required':False})
            atomic_state(pool_path,pool)
    def command(d):
        write_json(d/'actual_parameters.json',plan['parameters'])
        write_json(d/'run_plan.json',plan)
        shutil.copyfile(package/'provenance.json',d/'preflight_provenance.json')
        (d/'control').mkdir()
        frozen = d/'frozen_code'
        for name,digest in plan['code_sha256'].items():
            path=Path(name)
            if sha256_file(path)!=digest:raise ValueError('CODE_CHANGED_BEFORE_LAUNCH')
            target=frozen/path.relative_to(PROJECT_ROOT);target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(path,target)
        shutil.copyfile(PROJECT_ROOT/'vendor/Mirheo/LICENSE',d/'MIRHEO_LICENSE.txt')
        launch = ('#!/usr/bin/env bash\nset -euo pipefail\nsource '+shlex.quote(str(PROJECT_ROOT/'scripts/activate_mirheo.sh'))+
                  '\nexport PYTHONPATH='+shlex.quote(str(frozen))+
                  '\nexec /usr/bin/mpirun.openmpi --bind-to none -np 2 '+shlex.quote(source['environment']['python'])+
                  ' -B -u -m py_scripts.sdpd_diagnostics.equilibration_worker --spec '+shlex.quote(str(d/'actual_parameters.json'))+'\n')
        (d/'launch.sh').write_text(launch)
        return ['/bin/bash',str(d/'launch.sh')]
    d, record, hit = run_attempt(campaign,cc,c['case_id'],{'plan_sha256':plan['plan_sha256']},command,plan['allocation_s'],strict_allocation=True)
    if d is None:raise RuntimeError('RESOURCE_OR_BUDGET_BLOCKED_NO_LAUNCH '+str(record['status']))
    print('GPU_ATTEMPT '+json.dumps({'directory':str(d),'execution_status':record['status'],'cached':hit}),flush=True)
    return export(c,analyze=True)


def register_authorization(c, confirmation_file):
    """Explicit administrative action ONLY after real user approval; never called by preflight/execute."""
    confirmation_file = Path(confirmation_file).resolve()
    confirmation = read_json(confirmation_file)
    if confirmation.get('authority') != 'explicit_user_message':raise ValueError('EXPLICIT_USER_MESSAGE_REQUIRED')
    text_path = Path(confirmation['user_message_file']).resolve()
    message = text_path.read_text().strip()
    seconds = confirmation['additional_seconds']
    if not message or isinstance(seconds,bool) or not math.isfinite(seconds) or seconds<=0:
        raise ValueError('EXPLICIT_POSITIVE_SECONDS_AND_USER_MESSAGE_REQUIRED')
    package = export(c)
    request_path=package/'budget_request.json'
    request=read_json(request_path)
    if confirmation['budget_request_sha256']!=sha256_file(request_path) or confirmation['scope']!=request['authorization_scope']:
        raise ValueError('USER_CONFIRMATION_SCOPE_OR_REQUEST_MISMATCH')
    source=source_bundle(c)
    pool_path=Path(source['comparison_config']['shared_budget_pool'])
    record={'schema_version':1,'id':fingerprint(confirmation),'authority':confirmation['authority'],
            'user_message_file':str(text_path),'user_message_sha256':sha256_file(text_path),'user_message':message,
            'confirmation_file':str(confirmation_file),'confirmation_sha256':sha256_file(confirmation_file),
            'additional_seconds':seconds,'scope':confirmation['scope'],'single_task_limit_s':600,'gpu_concurrency':1,
            'budget_request_file':str(request_path),'budget_request_sha256':sha256_file(request_path),'recorded_at':now()}
    with exclusive_lock(read_json(pool_path)['gpu_lock']):
        pool=read_json(pool_path)
        from py_scripts.fluid_physics.budget_authorizations import authorization_records
        records=authorization_records(pool)
        if any(x['id']==record['id'] for x in records):raise ValueError('AUTHORIZATION_ALREADY_REGISTERED')
        d=output(Path(c['runs_root'])/c['campaign_id']/'authorizations'/record['id'])
        d.mkdir(parents=True,exist_ok=False)
        write_json(d/'shared_pool_before.json',pool)
        write_json(d/'authorization.json',record)
        pool.setdefault('authorization_records',[]).append({'path':str(d/'authorization.json'),'sha256':sha256_file(d/'authorization.json')})
        atomic_state(pool_path,pool)
    return d/'authorization.json'

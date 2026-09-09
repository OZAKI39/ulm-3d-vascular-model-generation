"""Extend the existing entry point, provenance, statistics and shared runner.

The current native SDPD restart contract is deficient. The funded three-process
test is diagnostic only; the formal segment loop remains fail-closed even if
particle averages happen to agree. No native changes or implicit authorization.
"""
import copy
import csv
import json
import math
from pathlib import Path
import shlex
import shutil
import numpy as np
import yaml
from py_scripts.fluid_physics.common import (PROJECT_ROOT,read_json,write_json,sha256_file,
    fingerprint,output,atomic_state,now)
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.runner import exclusive_lock,shared_budget_state,run_attempt
from .extended_restart import (native_contract,archive_checkpoint,verify_checkpoint,
                              compare_by_id,require_formal_restart,merge_segments)
from .extended_analysis import historical_supplement,observable_analysis,not_run


def json_safe(value):
    if isinstance(value,np.ndarray):return json_safe(value.tolist())
    if isinstance(value,np.generic):return json_safe(value.item())
    if isinstance(value,dict):return {str(k):json_safe(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [json_safe(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):raise ValueError('NONFINITE_JSON_VALUE')
    return value


def load_config(path):
    from .equilibration import load_config as legacy_load
    path=Path(path).resolve();raw=yaml.safe_load(path.read_text())
    if raw['schema_version']!=1 or raw['extra_authorized_gpu_seconds']!=0:
        raise ValueError('CONFIG_CANNOT_GRANT_GPU_AUTHORIZATION')
    base=legacy_load(output(raw['base_config']));c={**base,**raw}
    for key in ['data_root','runs_root','review_root','history_delivery_record','cpu_validation_record']:
        c[key]=str(output(c[key]))
    for key in ['campaign_id','case_id']:
        if not c[key] or '/' in c[key] or c[key] in ['.','..']:raise ValueError('INVALID_RUN_ID')
    p=c['plan'];r=c['restart_test']
    if (p['steps']!=800000 or p['dt_star']!=1e-6 or p['duration_star']!=.8
            or p['formal_window_star']!=[.6,.8] or p['diagnostic_window_star']!=[.4,.6]
            or p['sample_interval_star']!=.0002 or p['segment_steps']!=200000
            or p['segment_allocation_s']!=340 or p['checkpoint_mode']!='PingPong'
            or r['N1']!=2000 or r['N2']!=2000 or r['allocation_per_process_s']!=30 or r['process_count']!=3):
        raise ValueError('FIXED_EXTENDED_PLAN_CHANGED_REQUIRES_NEW_DESIGN')
    if c['policy']!=base['policy']:raise ValueError('FROZEN_SCOPE_POLICY_CHANGED')
    c['_config_path']=str(path);c['_config_sha256']=sha256_file(path)
    return c


def history(c):
    from .equilibration import source_bundle,verify_manifest,read_raw
    delivery_path=Path(c['history_delivery_record'])
    if sha256_file(delivery_path)!=c['history_delivery_sha256']:raise ValueError('HISTORY_DELIVERY_CHANGED')
    delivery=read_json(delivery_path);package=Path(delivery['package']);raw=Path(delivery['raw_directory'])
    # The old delivery also lists then-current editable Python files. Those are
    # historical code evidence, not immutable source files after this extension.
    for name,digest in delivery['artifact_sha256'].items():
        p=Path(name)
        if p.is_relative_to(PROJECT_ROOT/'py_scripts') or (p.parent==PROJECT_ROOT/'test_code'):continue
        if sha256_file(p)!=digest:raise ValueError('HISTORICAL_DELIVERY_ARTIFACT_CHANGED '+name)
    verify_package(package);verify_manifest(raw)
    plan=read_json(raw/'run_plan.json');spec=read_json(raw/'actual_parameters.json')
    if spec!=plan['parameters']:raise ValueError('HISTORICAL_FROZEN_PARAMETERS_CHANGED')
    if plan['criteria']!=c['criteria']:raise ValueError('OLD_CRITERIA_MUST_NOT_CHANGE')
    source=source_bundle(c)  # existing unit/native/library identity verification
    if spec['locked_units']!=source['spec']['locked_units'] or spec['candidate']!=source['spec']['candidate']:
        raise ValueError('FROZEN_PHYSICS_OR_UNITS_CHANGED')
    checks=read_json(package/'stationarity_checks.json')
    return {'delivery':delivery,'package':package,'raw':raw,'plan':plan,'spec':spec,
            'structure':checks['structure'],'data':read_raw(raw),'source':source,
            'execution':read_json(raw/'execution.json'),
            'frames':read_json(raw/'equilibration_observer.json')['snapshots']}


def build_plan(c,h,contract):
    from .equilibration import code_identity
    spec=copy.deepcopy(h['spec']);p=copy.deepcopy(c['plan']);dt=spec['dt_star']
    if dt!=p['dt_star']:raise ValueError('FROZEN_TIME_STEP_CHANGED')
    checkpoints=[str(x.relative_to(h['raw'])) for x in h['raw'].rglob('*') if x.is_file() and
                 (x.name.endswith('.state.txt') or x.name in ['state.txt','checkpoint_manifest.json'] or x.suffix in ['.xmf','.h5'])]
    if checkpoints:raise ValueError('HISTORICAL_CHECKPOINT_CANDIDATE_REQUIRES_EXPLICIT_PROVENANCE_AUDIT')
    spec.update(steps=p['steps'],desired_steps=p['steps'],desired_time_star=p['duration_star'],
                formal_window_star=p['formal_window_star'],source_category='FIXED_UNFORCED_EXTENDED_COLD_START',
                script_sha256=sha256_file(Path(__file__).with_name('extended_worker.py')),
                snapshot_marks_steps=[round(x/dt) for x in p['snapshot_times_star']])
    spec['task'].update(id=c['case_id'],desired_steps=p['steps'],desired_time_star=p['duration_star'],allocation_s=None)
    oldplan=h['plan'];rate=h['execution']['elapsed_monotonic_s']/oldplan['steps']
    tasks=[{'id':name,'kind':'restart_diagnostic','allocation_s':30,'steps':steps,'mode':mode} for name,steps,mode in
           [('restart_A',4000,'continuous'),('restart_B_save',2001,'save'),('restart_B_restore',2000,'diagnostic_restore')]]
    segments=[]
    for i in range(4):
        start=i*200000;end=(i+1)*200000
        segments.append({'segment_index':i,'start_step':start,'end_step':end,'absolute_time_star':end*dt,
                         'segment_steps':end-start,'worker_advance_steps':end-start+1,
                         'extra_checkpoint_trigger_steps':1,'allocation_s':340,
                         'parent_checkpoint':'NEW_COLD_INITIALIZATION' if not i else f'segment_{i-1:02d}_verified_checkpoint',
                         'actual_checkpoint_time_source':'native coordinator state file; reject if different from planned boundary'})
        tasks.append({'id':f'segment_{i:02d}','kind':'conditional_formal_segment','allocation_s':340,'steps':200001,'mode':'fresh' if i==0 else 'restore'})
    minimum=oldplan['criteria']['block_relaxation_multiples']*spec['domain_star'][1]**2/(4*math.pi**2*spec['candidate']['nu_prior_star'])
    plan={**p,'parameters':spec,'criteria':copy.deepcopy(c['criteria']),'restart_test':copy.deepcopy(c['restart_test']),
        'initialization':'COLD_START_NEW_TRAJECTORY','old_checkpoint_candidates':checkpoints,
        'old_snapshots_are_restart':False,'legal_increment_from_old_checkpoint_steps':None,
        'cold_start_committed_steps':800000,'cold_start_worker_steps_including_checkpoint_triggers':800004,
        'old_actual_time_star':.4,'old_actual_physical_time_s':.4*spec['locked_units']['t0'],
        'physical_duration_s':.8*spec['locked_units']['t0'],
        'minimum_block_duration_star':minimum,'minimum_block_samples':math.ceil(minimum/.0002),
        'planned_formal_samples':1000,'planned_samples_including_t0':4001,
        'segments':segments,'tasks':tasks,'scope_task_ids':[x['id'] for x in tasks],
        'cost_evidence':{'old_steps':oldplan['steps'],'old_elapsed_s':h['execution']['elapsed_monotonic_s'],
            'seconds_per_step':rate,'cold_start_proportional_s':800000*rate,
            'incremental_400000_s_if_valid_old_checkpoint_existed':400000*rate,
            'incremental_route_available':False,'unmeasured_storage_restart_overhead':True,
            'checkpoint_writes_per_200000_step_segment':1001,
            'checkpoint_storage_policy':'Native PingPong scratch limits disk growth; archive all referenced files after process exit into immutable segment versions. Existing 200-step u.run chunks reset the native checkpoint scheduler, causing a save at every chunk start.',
            'segment_expected_s_with_margin':200001*rate*p['safety_factor']+p['per_segment_setup_save_restore_output_reserve_s'],
            'restart_test_proportional_s':8001*rate,'restart_test_reserved_s':90},
        'maximum_total_allocation_s':1450,'formal_run_readiness':'BLOCKED_NATIVE_RNG_CHECKPOINT_CONTRACT',
        'budget_does_not_repair_native_contract':True,'restart_contract_sha256':fingerprint(contract),
        'native_environment':h['source']['environment'],
        'history_delivery_sha256':c['history_delivery_sha256'],'config_sha256':c['_config_sha256'],
        'code_sha256':code_identity(),'selection':None,'equilibrated_state_ready':False,
        'stop_rules':['plan end at 0.80*','restart state/RNG not validated','hash/ID/time mismatch','missing or duplicate samples',
                      'nonfinite state or illegal density','shared cumulative budget','task allocation and 600 s cap','no implicit retry or window relocation'],
        'single_long_process_fallback':{'enabled':False,'reason':'Would require a different explicitly approved single-task ceiling and design; it also cannot provide reliable SDPD RNG checkpoints in this build.'}}
    plan['plan_sha256']=fingerprint(plan)
    return plan


def budget_state(c,h,task_id='restart_A'):
    cc=copy.deepcopy(h['source']['comparison_config']);cc['read_only_unregistered_scope']=True
    return shared_budget_state(Path(c['runs_root'])/c['campaign_id'],cc,task_id=task_id)


def budget_request(c,plan,state):
    eligible=state['remaining_s'];additional=max(0.,1450-eligible)
    return {'status':'READY_TO_RUN_AWAITING_AUTHORIZATION' if additional else 'WITHIN_RECORDED_AUTHORIZATION',
        'authorization_scope':{'campaign_directory':str(Path(c['runs_root'])/c['campaign_id']),
                               'task_ids':plan['scope_task_ids'],'plan_sha256':plan['plan_sha256']},
        'total_requested_allocation_s':1450,'restart_comparison_allocation_s':90,
        'conditional_cold_start_allocation_s':1360,'per_segment_allocation_s':340,
        'expected_cold_start_s_with_margins':4*plan['cost_evidence']['segment_expected_s_with_margin'],
        'expected_total_upper_planning_s':4*plan['cost_evidence']['segment_expected_s_with_margin']+90,
        'current_scope_remaining_s':eligible,'original_remaining_s':state['original_remaining_s'],
        'old_scope_only_extension_remaining_s':sum(x['remaining_s'] for x in state['extension_accounts'] if x['campaign_directory']!=str(Path(c['runs_root'])/c['campaign_id'])),
        'additional_required_s':additional,'additional_request_whole_seconds':math.ceil(additional),
        'extra_authorized_gpu_seconds':0,'new_scope_recorded_extra_s':sum(r['additional_seconds'] for r in state['authorization_records'] if r['scope']=={'campaign_directory':str(Path(c['runs_root'])/c['campaign_id']),'task_ids':plan['scope_task_ids'],'plan_sha256':plan['plan_sha256']}),
        'single_task_limit_s':600,'GPU_concurrency':1,'user_approval_received':False,
        'maximum_time_star':.8,'maximum_physical_time_s':plan['physical_duration_s'],
        'scope_note':'One conditional request: three separate-process restore diagnostic tasks (30 s each), then four cold-start segments (340 s each) ONLY after independently verified native RNG persistence and real GPU restart PASS. Current native contract blocks all four formal segments; funding alone does not unblock them. Unspent balance remains in the same scoped account.',
        'automatic_extension':False,'budget_snapshot':state}


def require_budget(plan,state,task_ids=None):
    chosen=plan['scope_task_ids'] if task_ids is None else task_ids
    if not set(chosen)<=set(plan['scope_task_ids']):raise ValueError('TASK_OUTSIDE_APPROVED_PLAN')
    allocations=[t['allocation_s'] for t in plan['tasks'] if t['id'] in chosen]
    if any(a<=0 or a>600 for a in allocations):raise ValueError('SINGLE_TASK_LIMIT_EXCEEDED')
    if any(m['running_reservation_s'] for m in state['members']):raise RuntimeError('GPU_CONCURRENCY_RESERVATION_EXISTS')
    if state['remaining_s']+1e-9<sum(allocations):
        raise PermissionError('READY_TO_RUN_AWAITING_AUTHORIZATION: NO_GPU_STARTED; whole planned scope must fit before launch.')


def validation_record(c):
    from .equilibration import cpu_validation
    v=cpu_validation(c)
    if v and v.get('extended_config_sha256')!=c['_config_sha256']:raise ValueError('CPU_CONFIG_VALIDATION_STALE')
    return v


def export(c,analyze=False):
    h=history(c);contract=native_contract();plan=build_plan(c,h,contract);state=budget_state(c,h)
    request=budget_request(c,plan,state);v=validation_record(c)
    campaign=Path(c['runs_root'])/c['campaign_id']
    restart={'status':'RESTART_NOT_VALIDATED','execution_status':'NOT_RUN','evidence':None,
             'source_contract_blocker':contract['static_finding'],'formal_chain_allowed':False,'comparisons':None}
    if (campaign/'restart_validation.json').exists():
        restart=read_json(campaign/'restart_validation.json')
        if restart['plan_sha256']!=plan['plan_sha256']:raise ValueError('EXISTING_RESTART_TEST_DIFFERENT_PLAN_NO_RETRY')
    identity={'plan':plan['plan_sha256'],'request':fingerprint(request),'cpu':fingerprint(v) if v else None,
              'restart':fingerprint(restart),'raw':sha256_file(h['raw']/'output_sha256.json')}
    package=output(Path(c['data_root'])/('preparation_'+fingerprint(identity)[:16]))
    if package.exists():verify_package(package);return package
    supplement=historical_supplement(h['data'],h['plan'],h['structure'],h['raw'],h['frames'],h['execution']['elapsed_monotonic_s'])
    historical_status=observable_analysis(h['data'],h['plan'],h['structure'])
    statuses=not_run();statuses['RESTART_VALIDITY']={'status':'NOT_TESTED' if restart['execution_status']=='NOT_RUN' else 'FAIL','detail':restart}
    summary={'status':'READY_TO_RUN_AWAITING_AUTHORIZATION' if restart['execution_status']=='NOT_RUN' else 'RESTART_NOT_VALIDATED',
        'execution_status':'NOT_RUN','actual_steps':None,'actual_time_star':None,'actual_physical_time_s':None,
        'new_GPU_measurements':None,'initialization':'COLD_START_NEW_TRAJECTORY',
        'last_actual_historical_steps':400000,'last_actual_historical_time_star':.4,
        'last_actual_historical_physical_time_s':plan['old_actual_physical_time_s'],
        'restart_validity':restart['status'],'formal_run_readiness':plan['formal_run_readiness'],
        'CPU_validation':'PASS' if v else 'PENDING','production_selection':None,'selection':None,
        'equilibrated_state_ready':False,'human_review':'PENDING',
        'reason':'CPU preparation and historical supplemental analysis only. Nonzero-dt restore comparison awaits budget; formal continuation additionally lacks native SDPD RNG persistence.'}
    package.mkdir(parents=True)
    artifacts={'run_plan':plan,'budget_request':request,'restart_contract':contract,'restart_validation':restart,
        'segment_manifest':{'status':'NOT_RUN','actual_segments':[],'planned_segments':plan['segments'],'hashes':None,'absolute_time_star':None},
        'observable_status':statuses,'historical_observable_status':historical_status,
        'historical_supplement':supplement,'equilibration_summary':summary,
        'provenance':{'created_at':now(),'identity':identity,'history_delivery_record':c['history_delivery_record'],
            'history_delivery_sha256':c['history_delivery_sha256'],'history_raw':str(h['raw']),
            'history_raw_manifest_sha256':sha256_file(h['raw']/'output_sha256.json'),
            'native_environment':h['source']['environment'],'CPU_validation':v,'selection':None}}
    for name,value in artifacts.items():write_json(package/(name+'.json'),json_safe(value))
    shutil.copyfile(h['raw']/'raw_statistics.csv',package/'historical_raw_statistics.csv')
    from .equilibration_worker import RAW_COLUMNS
    with (package/'raw_statistics.csv').open('x') as f:csv.writer(f).writerow(RAW_COLUMNS)
    shutil.copyfile(c['_config_path'],package/'sdpd_equilibration_extended.yaml')
    (package/'report_zh.md').write_text(report_text(plan,request,supplement,historical_status,summary,package),encoding='utf-8')
    write_json(package/'package_sha256.json',{p.name:sha256_file(p) for p in package.iterdir() if p.is_file()})
    verify_package(package)
    preparation=output(campaign/'preparations'/package.name)
    if not preparation.exists():
        preparation.mkdir(parents=True)
        for name in ['run_plan.json','budget_request.json','restart_contract.json']:
            shutil.copyfile(package/name,preparation/name)
        write_json(preparation/'preparation_record.json',{'package':str(package),'manifest_sha256':sha256_file(package/'package_sha256.json'),'GPU_started':False})
    return package


def report_text(plan,b,s,old,summary,package):
    temp=s['details']['kBT_COM_star']['fixed_formal_description'];pressure=s['pressure_sampling'];cost=plan['cost_evidence']
    lines=['# 固定参数 SDPD：历史补充分析与恢复准备','',
        f"状态：{summary['status']}。本轮新 GPU 测量 NOT_RUN；CPU 测试 {summary['CPU_validation']}。",
        '长轨迹另有原生能力阻塞：RESTART_NOT_VALIDATED。预算获批也不能把当前粒子 checkpoint 当作完整随机过程恢复。selection=null；人工验收 PENDING。','',
        '## 来源与原参数','',
        f"配置固定引用已交付记录；逐一验证旧包、原始文件、报告和原生库哈希。旧实验 400000 步，t*=0.4，{plan['old_actual_physical_time_s']*1e6:.12f} µs；原结论 STATIONARITY_OR_SAMPLING_INCONCLUSIVE 保留。",
        '旧目录没有 state/XMF/HDF5 checkpoint。诊断快照不是续跑状态；规划独立冷启动，不拼接历史轨迹。',
        f"精确冻结参数：{json.dumps(plan['parameters']['candidate'],ensure_ascii=False)}",
        f"密度 1056 kg/m³；动力黏度 0.00345312 Pa·s；运动黏度 3.27e-6 m²/s；298.15 K；dt*={plan['dt_star']!r}，N=4096，盒 [8,8,8]。",
        f"单位：{json.dumps(plan['parameters']['locked_units'],ensure_ascii=False)}",
        '观测只复制/归约通道，没有添加力、额外恒温器或速度缩放。','',
        '## 旧数据的 CPU 补充证据','',
        f"正式窗口仍为 (0.30,0.40]，500 个样本。温度描述均值 {temp['mean']:.12g}*；四分窗均值 {temp['quarter_means']}；两半变化 0.789634%，超过原 0.5% 门槛。这不是可信稳态均值。",
        s['temperature_interpretation'],
        '温度：四段均值均下降，相关噪声修正下仍支持末段冷却；不能据此确定稳定偏差。HAC 的线性趋势、残差弱相关与渐近区间假设写入 JSON，未替换旧验收规则。',
        f"压力：每 200 步一个瞬时值，并非 200 步平均。原生应力包含保守、耗散及随机力。窗口标准差 {pressure['formal_sigma_pa']:.6f} Pa；门槛 {pressure['thermal_pressure_tolerance_pa']:.12f} Pa。",
        '压力四分窗均值不单调，多个 HAC 尺度的斜率区间均跨零：尚未证明存在持续漂移。旧同相位应力归约与 CPU 核密度核查已通过；没有删随机应力或用输入 EOS 代替机械压力。',
        f"即使独立同分布，500 个样本的平均压力 95% 半宽也约 {pressure['optimistic_iid_mean_halfwidth_pa']:.3f} Pa；这只是乐观算例，不是正式 CI。条件成本：{json.dumps(pressure['conditional_costs'],ensure_ascii=False)}",
        f"六个保存帧额外重建保守/耗散 virial，随机残差另列；与旧 CPU 重建差异最大 {s['pressure_components']['crosscheck_against_previous_CPU_density_reconstruction']['max_absolute_difference_star']:.6g}*。三个末段帧的保守压力升高约 {s['pressure_components']['late_conservative_virial_change_pa']:.6f} Pa，支持结构仍在松弛、背景被噪声遮盖的解释；三帧不能证明全窗口总压力漂移或给方差归因百分比。",
        '结构：核密度均值/方差持续下降，末段最近邻 q10 依次约 0.34867、0.36336、0.37626；范围 0.02759 rc 超过原 0.02 rc。三帧支持观察到分布变化，但不能给出充分独立的分布稳定置信度；方差缩小不等于结晶。',
        '独立信息：物理最小块长 40 个样本，500 点最多容纳 12 个完整块；这不是 12 个已验证独立块。过渡态 ACF 不被当作稳态 ACF，稳态 CI 留空。更晚 1000 点最多容纳 25 块，不能保证通过。',
        'PROPOSED（未应用）：评估逐步累计完整机械压力、仍每 200 步输出时间平均的开销和相关性；旧瞬时定义、压力门槛、旧结论均保留。','',
        '## 各指标分开判断（历史固定窗口）','',
        '| 指标 | 状态 | 说明 |','|---|---|---|']
    for key,row in old.items():lines.append(f"| {key} | {row['status']} | {row.get('reason','未测试；见记录')} |")
    lines+=['','### 保存帧机械压力分解（Pa，仍属于历史补充）','',
            '| 保存步数 | 保守 virial | 耗散 virial | 随机项及重建残差 | 原生总 virial |',
            '|---|---|---|---|---|']
    pressure_scale=plan['parameters']['locked_units']['si_per_star']['pressure']
    for f in s['pressure_components']['frames']:
        values=[f[k]*pressure_scale for k in ['conservative_virial_star','dissipative_virial_star','stochastic_plus_reconstruction_residual_star','total_native_virial_star']]
        lines.append('| '+str(f['step'])+' | '+' | '.join(f'{x:.6f}' for x in values)+' |')
    lines+=['','此表为 virial 分量；总机械压力还要加一次 COM 扣除后的动能项。没有删除任何随机贡献。',
        '新实验全部未测量；历史补充趋势不升级为新实验验收。','',
        '## 保存/恢复契约与实际限制','',
        '原生 SDPD 类及其基类没有 checkpoint/restart override，最终继承 MirObject 空实现；内核虽定义 RNG 序列化，但没有相互作用层调用。启用应力时还有复制内核的 RNG 与 StressManager 时间状态。DPD 的实现不能证明 SDPD 可恢复。当前构建的符号检查也只看到 DPD 与 MirObject 的相关实现。',
        '官方 interactions.py 用 DPD、dt=0；particle_vector.py 是 PV-only、dt=0；均不能代替目标非零 dt 测试。',
        '已实现退出后版本归档、物化软/硬链接、XMF 内部引用与 SHA 校验、完整标记、按 ID 比较、协调器步数与时间验证。文件完整性与随机状态完整性严格分开；缺时间、RNG、文件或哈希错误均拒绝正式恢复。',
        '每次 u.run 重建调度器并将 nExecutions 置零，所以启用 checkpoint 后每个 200 步块起点都会保存，不能把 checkpoint_every=200000 误当成全局周期。原生临时目录用 PingPong 限制体积，退出后物化成不可覆盖版本。每 20 万步段预计 1001 次保存，I/O 开销未实测，18 秒是预算余量而非测量承诺。边界另推进 1 步触发保存；归档以 state.txt 为准，额外步数计费且不混入正式样本。恢复不重抽速度、不扣 COM、不重置时间。',
        '已准备真实对照：A 连续 4000 步；B 保存 2000 步状态后退出（触发保存多推进的 1 步单列），另一个进程恢复并推进 2000 步。比较恢复前、首步及短段的粒子、时间、温度、压力、密度和动量。未运行不称 PASS；源代码 RNG 契约不通过时，即使均值接近也停止。',
        'CPU 测试只验证接口和保护逻辑，不替代真实 GPU 恢复。未修改或重编译原生库。','',
        '## 冻结后续设计和完整预算请求','',
        '正式窗口 (0.60,0.80]；趋势窗口 (0.40,0.60]；每 0.0002* 采样；最长 0.80*，不自动续到 1.0，也不移动窗口。冷启动 800000 个有效步，4 个 200000 步段，另外 4 步用于边界保存触发，均计费。',
        f"旧 400000 步实收 {cost['old_elapsed_s']:.12f} s；冷启动等比例 {cost['cold_start_proportional_s']:.6f} s；每段加 15% 和 18 s 未实测开销余量，估计 {cost['segment_expected_s_with_margin']:.6f} s，每段保留 340 s。",
        f"恢复对照 3×30=90 s；条件长实验 4×340=1360 s；一次请求总上限 1450 s。全局原授权 4093 s，已收 {b['budget_snapshot']['total_charged_or_reserved_s']:.12f} s；可用于此范围的旧余额 {b['current_scope_remaining_s']:.12f} s。旧 493 s 扩展已耗尽。",
        f"精确缺口 {b['additional_required_s']:.12f} s；申请整数追加 **{b['additional_request_whole_seconds']} GPU 秒**。并发 1，单任务不超过 600 s，实收可能更少，失败仍扣实际用量。",
        '该请求是有停止条件的完整上限，不是保证耗尽全部额度：当前原生 RNG 缺口意味着正式 1360 s 不可支出。默认只准备恢复诊断，失败或无法证明 RNG 连续就停止。预算不能单独解除原生能力阻塞；也不自动使用单个超 600 s 的长进程替代。',
        f"最大物理范围：{plan['physical_duration_s']*1e6:.12f} µs；所有新参数与生产 selection 仍为空。",
        '授权通过独立 --register-authorization 命令，绑定真实用户消息、此请求文件哈希与整组 task_ids。预检、查看、执行拒绝分支均不写授权、旧账本或启动 CUDA。','',
        '## 可执行命令','',
        '```bash','cd /home/lzy/projects/mirheo_starter',
        '.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --preflight-only',
        '.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --analyze-only',
        '.venv/bin/python -B -m test_code.review_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --open',
        '# 仅在真实批准并登记本请求后；当前代码会在恢复诊断后按原生缺口停止',
        '.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --execute','```','',
        '唯一下一建议：先完成有预算保护的真实非零 dt 恢复对照，确认此原生构建的恢复边界；恢复未验证前不启动新的长观察。']
    return '\n'.join(lines)+'\n'


def launch_task(c,h,plan,task,spec,cc):
    """One shared runner call; the same campaign account covers every segment."""
    campaign=Path(c['runs_root'])/c['campaign_id']
    def command(d):
        write_json(d/'actual_parameters.json',spec);write_json(d/'run_plan.json',plan);(d/'control').mkdir()
        frozen=d/'frozen_code'
        for name,digest in plan['code_sha256'].items():
            src=Path(name)
            if sha256_file(src)!=digest:raise ValueError('CODE_CHANGED_BEFORE_GPU_LAUNCH')
            dst=frozen/src.relative_to(PROJECT_ROOT);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        shutil.copyfile(PROJECT_ROOT/'vendor/Mirheo/LICENSE',d/'MIRHEO_LICENSE.txt')
        launch='#!/usr/bin/env bash\nset -euo pipefail\nsource '+shlex.quote(str(PROJECT_ROOT/'scripts/activate_mirheo.sh'))+'\nexport PYTHONPATH='+shlex.quote(str(frozen))+'\nexec /usr/bin/mpirun.openmpi --bind-to none -np 2 '+shlex.quote(h['source']['environment']['python'])+' -B -u -m py_scripts.sdpd_diagnostics.extended_worker --spec '+shlex.quote(str(d/'actual_parameters.json'))+'\n'
        (d/'launch.sh').write_text(launch);return ['/bin/bash',str(d/'launch.sh')]
    return run_attempt(campaign,cc,task['id'],{'plan':plan['plan_sha256'],'task':task['id']},command,task['allocation_s'],strict_allocation=True)


def execute(c):
    from .equilibration import cpu_validation
    package=export(c);plan=read_json(package/'run_plan.json');request=read_json(package/'budget_request.json')
    h=history(c);campaign=Path(c['runs_root'])/c['campaign_id'];contract=read_json(package/'restart_contract.json')
    if (campaign/'restart_validation.json').exists():return package
    require_budget(plan,budget_state(c,h))  # before any runner mutation/CUDA/resource probe
    cpu_validation(c,required=True);validation_record(c)
    cc=copy.deepcopy(h['source']['comparison_config']);cc['campaign_id']=c['campaign_id']
    pool_path=Path(cc['shared_budget_pool'])
    with exclusive_lock(read_json(pool_path)['gpu_lock']):
        state=budget_state(c,h);require_budget(plan,state)
        if not any(r['scope']==request['authorization_scope'] for r in state['authorization_records']):
            raise PermissionError('EXACT_SCOPED_AUTHORIZATION_REQUIRED_NO_GPU_STARTED')
        pool=state['pool']
        if not any(Path(m['directory']).resolve()==campaign.resolve() for m in pool['members']):
            write_json(campaign/'shared_pool_before_registration.json',pool)
            pool['members'].append({'directory':str(campaign),'campaign_id':c['campaign_id'],'ledger_required':False})
            atomic_state(pool_path,pool)
    outputs={};parent=None
    for task in plan['tasks'][:3]:
        spec=copy.deepcopy(plan['parameters'])
        spec.update(steps=task['steps'],extended_task=task,extended_plan=plan,
                    restore_directory=str(parent) if task['mode']=='diagnostic_restore' else None,
                    restore_manifest_sha256=sha256_file(parent/'checkpoint_manifest.json') if task['mode']=='diagnostic_restore' else None)
        d,record,_=launch_task(c,h,plan,task,spec,cc)
        if d is None or record['status']!='COMPLETED':
            write_json(campaign/'restart_validation.json',{'status':'RESTART_NOT_VALIDATED','execution_status':'FAILED_OR_PARTIAL',
                'plan_sha256':plan['plan_sha256'],'evidence':'REAL_GPU_SEPARATE_PROCESSES','failed_task':task['id'],'formal_chain_allowed':False})
            return export(c,analyze=True)
        outputs[task['id']]=str(d)
        if task['mode']=='save':
            parent=archive_checkpoint(d/'native_checkpoint',campaign/'checkpoints'/'restart_B',spec=plan['parameters'],
                contract=contract,completion=read_json(d/'worker_completion.json'))
    validation=validate_real_restart(outputs,plan,contract)
    write_json(campaign/'restart_validation.json',json_safe(validation))
    # This source version deliberately has no successful route through this gate.
    # A new native contract would need separate user-authorized native work.
    if validation['status']=='PASS':
        require_formal_restart(contract,validation)
        run_formal_segments(c,plan,h,contract,validation,cc,package)
    return export(c,analyze=True)


def validate_real_restart(outputs,plan,contract):
    paths={k:Path(v) for k,v in outputs.items()};r=plan['restart_test'];N1=r['N1'];N2=r['N2'];checks={}
    from .equilibration import verify_manifest
    for directory in paths.values():verify_manifest(directory)
    def state(task,name):
        with np.load(paths[task]/(name+'.npz'),allow_pickle=False) as z:return {k:z[k] for k in z.files}
    a=state('restart_A','initial_state');b=state('restart_B_save','initial_state')
    checks['identical_fresh_initial_state']=compare_by_id(a,b,plan['parameters']['domain_star'],
        position_atol=r['preadvance_position_atol_star'],velocity_atol=r['preadvance_velocity_atol_star'])
    checks['identical_fresh_initial_state']['time_correct']=int(a['step'])==int(b['step'])==0 and float(a['time_star'])==float(b['time_star'])==0
    metadata=[read_json(p/'initialization_record.json') for p in paths.values()]
    checks['same_rank_layout_and_initialization']=all(m['rank_layout']==[1,1,1] and m['world_size']==2 for m in metadata) and [m['velocity_initialization_count'] for m in metadata]==[1,1,0] and [m['COM_subtraction_count'] for m in metadata]==[1,1,0]
    checks['same_physics_and_library']=all(read_json(p/'run_plan.json')==plan for p in paths.values())
    for name,step,pa,va in [('preadvance',N1,r['preadvance_position_atol_star'],r['preadvance_velocity_atol_star']),
                            ('first_step',N1+1,r['first_step_position_atol_star'],r['first_step_velocity_atol_star']),
                            ('short',N1+N2,r['short_position_atol_star'],r['short_velocity_atol_star'])]:
        a=state('restart_A',f'state_{step}');b=state('restart_B_restore','restored_before_advance' if name=='preadvance' else f'state_{step}')
        checks[name]=compare_by_id(a,b,plan['parameters']['domain_star'],position_atol=pa,velocity_atol=va)
        checks[name]['time_correct']=int(a['step'])==int(b['step'])==step and abs(float(a['time_star'])-float(b['time_star']))<=2e-10
        checks[name]['mass_correct']=float(a['mass_star'])==float(b['mass_star'])==plan['parameters']['m_star']
        if name!='preadvance':
            aa=read_json(paths['restart_A']/f'observation_{step}.json');bb=read_json(paths['restart_B_restore']/f'observation_{step}.json')
            limits={'kBT_COM_star':r['temperature_jump_atol_star'],'pressure_star':r['pressure_jump_atol_star'],
                    'kernel_number_mean':r['kernel_mean_jump_atol_star'],**{'mean_v'+axis:r['COM_jump_atol_star'] for axis in 'xyz'}}
            checks[name]['observables']={k:{'absolute_difference':abs(aa[k]-bb[k]),'limit':tol,
                                           'status':'PASS' if abs(aa[k]-bb[k])<=tol else 'FAIL'} for k,tol in limits.items()}
    from .equilibration import read_raw
    A=read_raw(paths['restart_A']);B=read_raw(paths['restart_B_save']);C=read_raw(paths['restart_B_restore'])
    sequence=np.concatenate([B['step'][B['step']<=N1],C['step']])
    checks['sample_sequence_pass']=bool(np.array_equal(sequence,A['step']))
    numerical=all(checks[k]['status']=='PASS' and checks[k]['time_correct'] and checks[k].get('mass_correct',True)
        and all(x['status']=='PASS' for x in checks[k].get('observables',{}).values())
        for k in ['identical_fresh_initial_state','preadvance','first_step','short'])
    numerical=numerical and all(checks[k] for k in ['sample_sequence_pass','same_rank_layout_and_initialization','same_physics_and_library'])
    passed=numerical and contract['native_rng_persistence_proven']
    return {'status':'PASS' if passed else 'RESTART_NOT_VALIDATED','execution_status':'COMPLETED_COMPARISON',
        'evidence':'REAL_GPU_SEPARATE_PROCESSES','plan_sha256':plan['plan_sha256'],'numerical_comparison_pass':bool(numerical),
        'RNG_persistence_pass':contract['native_rng_persistence_proven'],'checks':checks,'formal_chain_allowed':bool(passed),
        'outputs':{k:{'directory':str(p),'manifest_sha256':sha256_file(p/'output_sha256.json')} for k,p in paths.items()}}


def run_formal_segments(c,plan,h,contract,validation,cc,package):
    """Bounded segment loop, unusable until native RNG persistence is demonstrated.

    The current build fails the first gate. The loop is real orchestration, not
    CPU simulation; it never treats a particle archive as a full checkpoint.
    """
    require_formal_restart(contract,validation)
    from .equilibration import verify_manifest
    campaign=Path(c['runs_root'])/c['campaign_id'];records=[];parent=None;parent_digest=None
    require_budget(plan,budget_state(c,h,task_id='segment_00'),[t['id'] for t in plan['tasks'][3:]])
    try:
        for segment,original_task in zip(plan['segments'],plan['tasks'][3:]):
            require_formal_restart(contract,validation)
            if parent:verify_checkpoint(parent,require_rng=True,expected_manifest=parent_digest)
            task={**original_task,'committed_end_step':segment['end_step']}
            spec=copy.deepcopy(plan['parameters']);spec.update(steps=task['steps'],extended_task=task,extended_plan=plan,
                restore_directory=str(parent) if parent else None,restore_manifest_sha256=parent_digest)
            d,execution,_=launch_task(c,h,plan,task,spec,cc)
            if d is None or execution['status']!='COMPLETED':raise RuntimeError('FORMAL_SEGMENT_INCOMPLETE_NO_NEXT_SEGMENT')
            verify_manifest(d);completion=read_json(d/'worker_completion.json')
            if completion['status']!='COMPLETED_PLANNED_STEPS' or completion['steps']!=segment['end_step']+1:
                raise ValueError('SEGMENT_ACTUAL_ADVANCE_MISMATCH')
            checkpoint=archive_checkpoint(d/'native_checkpoint',campaign/'checkpoints'/task['id'],
                spec=plan['parameters'],contract=contract,completion=completion,parent=parent_digest)
            cm=verify_checkpoint(checkpoint,require_rng=True)
            if cm['absolute_step']!=segment['end_step']:raise ValueError('SCHEDULED_CHECKPOINT_NOT_AT_REQUESTED_BOUNDARY')
            with (d/'raw_statistics.csv').open() as stream:
                rows=[{k:(float(v) if v else None) for k,v in r.items()} for r in csv.DictReader(stream)]
            for row in rows:row['step']=int(row['step'])
            digest=sha256_file(checkpoint/'checkpoint_manifest.json')
            record={'status':'COMPLETED',**segment,'directory':str(d),'parent_checkpoint':parent_digest,
                    'checkpoint_directory':str(checkpoint),'checkpoint_sha256':digest,
                    'output_manifest_sha256':sha256_file(d/'output_sha256.json'),'rows':rows,
                    'absolute_step':cm['absolute_step'],'absolute_time_star':cm['absolute_time_star'],
                    'actual_advanced_steps':completion['segment_advanced_steps'],'actual_charge_s':execution['elapsed_monotonic_s']}
            records.append(record)
            # Each prefix must form a complete chain up to its actual boundary.
            merge_segments(records,plan['dt_star'],plan['parameters']['snapshot_every'],segment['end_step'])
            write_json(campaign/f'segment_record_{segment["segment_index"]:02d}.json',json_safe(record))
            parent,parent_digest=checkpoint,digest
        rows=merge_segments(records,plan['dt_star'],plan['parameters']['snapshot_every'],plan['steps'])
        with (campaign/'merged_raw_statistics.csv').open('x') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        write_json(campaign/'segment_manifest.json',{'status':'COMPLETED','plan_sha256':plan['plan_sha256'],
            'actual_segments':records,'absolute_step':plan['steps'],'absolute_time_star':plan['duration_star'],
            'merged_raw_sha256':sha256_file(campaign/'merged_raw_statistics.csv')})
    except (ValueError,RuntimeError,OSError) as exc:
        write_json(campaign/'segment_manifest.json',json_safe({'status':'PARTIAL_OR_RESTART_INVALID',
            'plan_sha256':plan['plan_sha256'],'actual_segments':records,'stop_reason':str(exc),
            'absolute_step':records[-1]['absolute_step'] if records else 0,'formal_statistics':None}))
        raise

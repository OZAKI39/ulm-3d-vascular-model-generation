"""Read-only reconstruction of completed extended GPU diagnostics.

Never changes frozen worker code, authorization, raw data, or restart verdicts.
Analysis revisions have their own identity, separate from the executed plan.
"""
import copy
import csv
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT, read_json, sha256_file, fingerprint
from .extended_restart import verify_checkpoint, checkpoint_channel_forms, compare_by_id, inside


def executed_plan(campaign, proposed, terminal):
    from .equilibration import verify_manifest
    ledger = read_json(Path(campaign)/'budget_ledger.json')
    if not ledger['attempts']: raise ValueError('TERMINAL_RESTART_WITHOUT_ATTEMPTS')
    first = Path(ledger['attempts'][0]['directory'])
    if first.parent.resolve() != Path(campaign).resolve(): raise ValueError('ATTEMPT_OUTSIDE_CAMPAIGN')
    verify_manifest(first)
    plan = read_json(first/'run_plan.json')
    if fingerprint({k:v for k,v in plan.items() if k!='plan_sha256'}) != plan['plan_sha256']:
        raise ValueError('EXECUTED_PLAN_DIGEST_MISMATCH')
    if terminal['plan_sha256'] != plan['plan_sha256']: raise ValueError('TERMINAL_PLAN_MISMATCH')
    def design(p):
        p = copy.deepcopy(p)
        for key in ['plan_sha256','code_sha256']: p.pop(key, None)
        p['parameters'].pop('script_sha256', None)
        return p
    if design(plan) != design(proposed):
        raise ValueError('EXECUTED_DESIGN_CHANGED_REANALYSIS_ONLY')
    return plan


def execution_evidence(campaign, plan, terminal):
    from .equilibration import verify_manifest
    root = Path(campaign); ledger = read_json(root/'budget_ledger.json'); jobs = []
    tasks = {t['id']:t for t in plan['tasks']}; seen = set()
    for attempt in ledger['attempts']:
        task = attempt['task_id']; directory = Path(attempt['directory'])
        if task not in tasks or task in seen or tasks[task]['kind'] != 'restart_diagnostic':
            raise ValueError('UNEXPECTED_OR_RETRIED_DIAGNOSTIC_ATTEMPT')
        seen.add(task)
        if directory.parent.resolve()!=root.resolve() or attempt['cache_identity']!={'plan':plan['plan_sha256'],'task':task}:
            raise ValueError('ATTEMPT_PROVENANCE_MISMATCH')
        verify_manifest(directory)
        if read_json(directory/'run_plan.json') != plan: raise ValueError('ATTEMPT_PLAN_MISMATCH')
        raw = read_json(directory/'execution.json')
        if raw['status'] != attempt['status'] or raw['elapsed_monotonic_s'] != attempt['charged_s']:
            raise ValueError('ATTEMPT_LEDGER_EXECUTION_MISMATCH')
        complete = read_json(directory/'worker_completion.json') if (directory/'worker_completion.json').is_file() else None
        samples = []
        if (directory/'raw_statistics.csv').is_file():
            with (directory/'raw_statistics.csv').open() as stream:
                samples = [{k:float(v) if v else None for k,v in r.items()} for r in csv.DictReader(stream)]
        jobs.append({'task_id':task, 'directory':str(directory), 'status':attempt['status'],
            'charged_s':attempt['charged_s'], 'allocation_s':attempt['reserved_s'], 'exit_code':raw['exit_code'],
            'timeout':raw['timeout'], 'remaining_own_group_processes':raw['remaining_own_group_processes'],
            'actual_steps':complete['steps'] if complete else None,
            'actual_time_star':complete['time_star'] if complete else None,
            'actual_physical_time_s':complete['time_star']*plan['parameters']['locked_units']['t0'] if complete else None,
            'completion':complete, 'raw_rows':len(samples), 'samples':samples,
            'output_manifest_sha256':sha256_file(directory/'output_sha256.json'),
            'initialization_record':read_json(directory/'initialization_record.json') if (directory/'initialization_record.json').exists() else None,
            'restored_preadvance_observed':(directory/'restored_before_advance.npz').is_file()})
    if not jobs or any(j['status']=='RUNNING' for j in jobs): raise ValueError('TERMINAL_WITH_UNFINISHED_ATTEMPTS')
    return {'status':'RESTART_DIAGNOSTIC_FAILED' if terminal.get('failed_task') else 'RESTART_DIAGNOSTIC_NOT_VALIDATED',
        'evidence':'REAL_GPU_SEPARATE_PROCESSES', 'plan_sha256':plan['plan_sha256'],
        'jobs':jobs, 'GPU_attempt_count':len(jobs), 'charged_s':sum(j['charged_s'] for j in jobs),
        'formal_execution_status':'NOT_RUN', 'formal_segments':[], 'automatic_retry':False,
        'terminal_record':str(root/'restart_validation.json'),
        'terminal_record_sha256':sha256_file(root/'restart_validation.json'),
        'campaign_ledger_sha256':sha256_file(root/'budget_ledger.json'),
        'equilibrium_evidence':False, 'independent_trajectories_are_not_concatenated':True}


def checkpoint_failure_audit(campaign, plan, terminal):
    """Identify actual reader failure; compare saved HDF state on CPU by ID."""
    root = Path(campaign); directory = root/'checkpoints/restart_B'
    if not directory.exists():
        return {'status':'NO_COMPLETED_CHECKPOINT', 'restore_checks':'NOT_REACHED'}
    manifest = verify_checkpoint(directory, require_rng=False)
    if manifest['frozen_parameters_sha256'] != fingerprint(plan['parameters']):
        raise ValueError('ARCHIVE_PARAMETERS_DIFFER_FROM_EXECUTED_PLAN')
    formats = checkpoint_channel_forms(directory)
    import h5py  # Read-only CPU inspection; no import of Mirheo or CUDA.
    datasets = {}
    for name in manifest['files']:
        if name.endswith('.h5'):
            with h5py.File(directory/name, 'r') as f:
                datasets[name] = {key:{'shape':list(f[key].shape),'dtype':str(f[key].dtype)} for key in f}
    xmf = ET.parse(directory/'pv.PV.xmf')
    def dataset(node):
        ref, key = (node.text or '').strip().split(':',1)
        with h5py.File(inside(directory,ref), 'r') as stream: return stream[key][:]
    attrs = {a.get('Name'):a.find('DataItem') for a in xmf.iter('Attribute')}
    state = {'positions':dataset(xmf.find('.//Geometry/DataItem')),
             'velocities':dataset(attrs['velocities']), 'ids':dataset(attrs['ids']).reshape(-1)}
    def saved(task, name):
        with np.load(root/task/(name+'.npz'), allow_pickle=False) as z: return dict(z)
    r = plan['restart_test']; domain = plan['parameters']['domain_star']; n1 = r['N1']
    def compare(a,b):
        return compare_by_id(a,b,domain,position_atol=r['preadvance_position_atol_star'],
                             velocity_atol=r['preadvance_velocity_atol_star'])
    comparisons = {
        'identical_fresh_initial_state':compare(saved('restart_A','initial_state'),saved('restart_B_save','initial_state')),
        'checkpoint_HDF_vs_own_saved_boundary':compare(saved('restart_B_save',f'state_{n1}'),state),
        'independent_A_vs_B_before_any_restore':compare(saved('restart_A',f'state_{n1}'),saved('restart_B_save',f'state_{n1}'))}
    # The original coordinator comparison was A(N1) versus restored B(N1).
    # The observed pre-save A/B difference is already above the frozen velocity
    # tolerance. Do not silently relax it or call the entire restore comparison PASS.
    comparisons['independent_A_vs_B_before_any_restore']['interpretation'] = (
        'Pre-existing A/B numerical difference, before restore. Report against the frozen 2e-6 '
        'preadvance limits; not a measured restart-induced jump and not permission to relax limits.')
    failure_task = terminal.get('failed_task')
    console = root/failure_task/'console.log' if failure_task else None
    log = root/failure_task/'mirheo_00000.log' if failure_task else None
    evidence = []
    for path in [console,log]:
        if path and path.exists():
            evidence.extend({'path':str(path),'line':i,'text':line} for i,line in enumerate(path.read_text(errors='replace').splitlines(),1)
                            if 'Unrecognised form Other' in line)
    native = PROJECT_ROOT/'vendor/Mirheo'
    files = ['src/mirheo/core/xdmf/type_map.h','src/mirheo/core/xdmf/channel.cpp',
             'src/mirheo/core/xdmf/xmf_helpers.cpp','src/mirheo/plugins/particle_channel_saver.cpp',
             'src/mirheo/core/integrators/vv.cu']
    return {'status':'RESTART_NOT_VALIDATED', 'inspection':'CPU_READ_ONLY_AFTER_REAL_GPU_FAILURE',
        'archive_directory':str(directory), 'archive_manifest_sha256':sha256_file(directory/'checkpoint_manifest.json'),
        'file_integrity':'PASS', 'native_readability':formats, 'HDF_datasets':datasets,
        'checkpoint_absolute_step':manifest['absolute_step'], 'checkpoint_absolute_time_star':manifest['absolute_time_star'],
        'checkpoint_absolute_physical_time_s':manifest['absolute_time_star']*plan['parameters']['locked_units']['t0'],
        'committed_boundary_step':n1, 'save_worker_advanced_steps':n1+1,
        'boundary_vs_trigger_note':'Saved state is step 2000; step 2001 only triggers native save, and is separately charged.',
        'comparisons':comparisons,
        'restored_preadvance_state':'NOT_REACHED' if not (root/'restart_B_restore/restored_before_advance.npz').exists() else 'OBSERVED',
        'first_step_comparison':'NOT_RUN', 'short_segment_comparison':'NOT_RUN',
        'temperature_pressure_density_momentum_restore_jumps':'NOT_MEASURED',
        'restore_sample_sequence':'NOT_TESTED', 'RNG_persistence_pass':False,
        'failure_evidence':evidence,
        'direct_cause':'saved_forces has Typeinfo Other/Datatype Force; native readDataSet aborts on Other. HDF has one column, not three force components.',
        'Python_issue':'Diagnostic Force channel was marked persistent without checking native XMF serialization support; byte integrity alone incorrectly allowed native restore to begin.',
        'CPU_fix':'Reject unsupported persistent channel forms before Mirheo/CUDA initialization and before coordinator.restart. Preserve all original files; do not relabel, drop, or fake missing force data.',
        'not_fixed':'This guard does not implement valid Force serialization or native SDPD RNG persistence; no GPU retry was performed.',
        'independent_blocker':'SDPD interaction inherits no-op checkpoint/restart; kernel and stress-wrapper RNG state is not persisted. Fixing the XMF reader failure alone is insufficient.',
        'channel_phase':'Saved diagnostic channels belong to the preceding force step. Position/velocity/ID HDF match the B save boundary; ephemeral force/density must be recomputed on advance.',
        'native_source_evidence':[{'path':str(native/name),'sha256':sha256_file(native/name)} for name in files],
        'formal_chain_allowed':False, 'selection':None}


def terminal_report(plan, budget, summary, execution, audit, historical_report):
    comparisons = audit.get('comparisons',{}); before = comparisons.get('independent_A_vs_B_before_any_restore',{})
    lines = ['# 固定参数 SDPD：真实恢复测试失败与 CPU 原因核查', '',
        '**RESTART_NOT_VALIDATED。三进程恢复诊断已执行，恢复读取失败，长实验按约定停止。**',
        f"本轮计费 {execution['charged_s']:.12f} GPU 秒；CPU 修复验证 {summary['CPU_validation']}。没有自动重试。selection=null；人工验收 PENDING。", '',
        '## 实际执行与绝对时间', '',
        '| 独立任务 | 结果 | 实际步数 | 实际 t* | GPU 秒 |', '|---|---|---:|---:|---:|']
    for j in execution['jobs']:
        lines.append(f"| {j['task_id']} | {j['status']} (exit {j['exit_code']}) | {j['actual_steps']} | {j['actual_time_star']} | {j['charged_s']:.9f} |")
    lines += ['', 'A、B 分别从同一初态冷启动，不能拼成连续轨迹；此短对照仅检验恢复，不检验液体平衡。A 达到约 t*=0.004；B 保存的真实边界为第 2000 步、t*=0.002，随后多推进 1 步触发写盘。',
        'B-restore 在 u.restart 内退出，尚未得到恢复后推进前状态；不能把计划再推进的 2000 步记为已完成。三进程均未超时，runner 未发现遗留本任务进程。',
        '本轮正式长轨迹 NOT_RUN，(0.60,0.80] 没有新样本。最近完成的液体长观察仍是旧 400000 步、t*=0.40、12.662309354 µs，旧结论与门槛保留。', '',
        '## 直接失败原因及独立阻塞', '',
        '原生日志：`Unrecognised form Other`，定位到 `core/xdmf/xmf_helpers.cpp:103`。实际 `pv.PV.xmf` 中 `saved_forces` 为 `Other / Force`；HDF5 实际形状是 `[4096, 1]`。`type_map.h` 未给 Force 定义 DataForm，回退到 Other；`channel.cpp` 把它当成 1 分量写出，而读取器拒绝 Other。',
        '`saved_stresses` 是受支持的 Tensor6 / Stress，HDF 形状 `[4096, 6]`，不是这次错误源。将 saved_forces 标签改成 Vector 不能补回不存在的三维数据；本轮没有改写原始 checkpoint。',
        '项目 Python 的问题是：启用了无法原生序列化的持久化诊断 Force 通道，却仅检查文件哈希。已先用失败测试复现，再增加 CPU 通道格式检查，在导入 Mirheo/CUDA 和调用 restart 前拒绝该文件。该修复是明确拒绝无效输入，不代表 Force 序列化已修好。',
        '另一项独立限制仍然成立：当前 SDPD 相互作用没有接入内核随机状态的 checkpoint/restart，应力包装器还有独立 RNG 与调度状态。即使解决文件读取，也不能据此建立完整随机过程续跑。未修改、重编译或升级原生库。', '',
        '## 保存状态与恢复比较', '',
        '| 核查 | 状态 | 证据 |', '|---|---|---|',
        '| 文件完整性 | PASS | manifest、内部引用和每个文件 SHA 验证 |',
        '| 相同冷启动初态 | PASS | 按 int64 粒子 ID 比较，位置/速度最大差均为 0 |',
        '| HDF 与 B 自身第 2000 步 | PASS | 按 ID 比较，位置/速度最大差均为 0；使用原 2e-6 门槛 |',
        f"| 独立 A/B 第 2000 步 | {before.get('status')} | 最大速度差 {before.get('max_velocity_error_star')}*，原门槛 2e-6*；已在恢复前出现 |",
        '| 恢复后推进前、第一步、短段 | NOT_RUN / NOT_REACHED | 原生读取失败，缺实测状态 |',
        '| 温度/压力/密度/动量恢复跳变 | NOT_MEASURED | 不用均值接近或合成数据替代 |',
        '| SDPD RNG 连续性 | 未验证 | 原生持久化契约缺失 |', '',
        'A/B 的保存前速度差已超过原先用于 A 对恢复 B 的严格门槛；它不是恢复造成的跳变，也不能忽略。这里只补充描述，未事后放宽门槛。以后验证应单独区分“恢复状态对其自身保存状态”和“独立进程演化差异”，需要在新试验前冻结设计。', '',
        '## 授权、停止和余额', '',
        f"批准已登记：追加 {budget['extra_authorized_gpu_seconds']:g} 秒，绑定原计划 {plan['plan_sha256']} 与指定 7 个任务。原 3600 秒基额和旧 493 秒扩展保留，总授权 {budget['budget_snapshot']['limit_s']:.0f} 秒。",
        f"全局累计计费 {budget['budget_snapshot']['total_charged_or_reserved_s']:.12f} 秒；可用于本范围的余额 {budget['current_scope_remaining_s']:.12f} 秒，其中原基额余额 {budget['original_remaining_s']:.12f} 秒。",
        '原批准范围上限为 3×30 + 4×340 = 1450 秒；失败后不支出 1360 秒条件长段。未花余额留在原用途账户，不清空、不回填、也不自动授权重试。本次新增申请为 0 秒。',
        '并发 1，单任务上限 600 秒；没有改成单个超上限长进程。正式窗口仍为 (0.60,0.80]，趋势窗口仍为 (0.40,0.60]，最长 0.80*。', '',
        '## 可复现的只读分析与页面', '', '```bash', 'cd /home/lzy/projects/mirheo_starter',
        '.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --analyze-only',
        '.venv/bin/python -B -m test_code.review_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --open', '```',
        '已有 terminal 失败记录时，--execute 也只返回分析包，不重新启动 GPU。执行计划和各进程 frozen_code 保留原哈希；CPU 分析修订另记代码身份，不把修订后的代码冒充此次 GPU 所用版本。',
        '读取失败后自动生成的旧 preparation_779eb32f8aae47f0 包中，NOT_RUN 总状态和重复预算请求是 Python 汇报错误。原包保留，本报告及 restart_execution.json 对其作明确更正。', '',
        '下一优先项：先解决持久化通道与 SDPD 随机状态的完整恢复契约，再做新的非零 dt 验证。当前约束禁止原生修改，因此本轮在 CPU 证据和拒绝保护处停止，不能宣称可靠恢复或后段平衡已完成。', '',
        '## 历史物理结果的补充证据（不属于本次短诊断）', '']
    # Preserve the existing source/physics and historical sections, without its
    # superseded pending execution and budget prose.
    source = historical_report.split('## 来源与原参数',1)[1].split('## 保存/恢复契约与实际限制',1)[0]
    source = source.replace('新实验全部未测量；历史补充趋势不升级为新实验验收。',
                            '本轮正式液体长观察未测量；短恢复诊断不作为物性证据，历史补充趋势不升级为新实验验收。')
    return '\n'.join(lines)+'\n'+source

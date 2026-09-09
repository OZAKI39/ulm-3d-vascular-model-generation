"""Candidate-local decisions, two distinct cost comparisons and immutable reports."""
from pathlib import Path
import csv
import json
import math
import shutil
import yaml
from py_scripts.fluid_physics.common import PROJECT_ROOT,output,read_json,write_json,sha256_file,fingerprint,now,environment_identity
from py_scripts.fluid_physics.units import Units
from py_scripts.fluid_physics.runner import shared_budget_state
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.reporting import eos_analysis,sensitivity
from .models import frozen_targets,native_contract
from .experiments import verify_historical
from .analysis import analyze_run


def choose(assessments):
    qualified=[a for a in assessments if a['status']=='QUALIFIED_PROPOSED']
    return {'method':qualified[0]['method'],'candidate_id':qualified[0]['candidate_id']} if len(qualified)==1 else None


def assess(tasks,case,u,c):
    identities=sorted({(t['method'],t['candidate_id'],t['comparison_group']) for t in tasks})
    rows=[];tol=c['proposed_acceptance']
    for method,cid,group in identities:
        ts=[t for t in tasks if (t['method'],t['candidate_id'],t['comparison_group'])==(method,cid,group)]
        flows=[t for t in ts if t['test_kind'] in ('flow','half_dt','half_force')]
        equilibrium=[t for t in ts if t['test_kind']=='equilibrium']
        base=next((t for t in flows if t['test_kind']=='flow'),None)
        eos=eos_analysis(ts,case,u,tol);sens=sensitivity(ts,tol)
        gates={};reasons=[]
        gates['interface']=any(t['execution'].get('exit_code')==0 and t.get('pressure_channels_verified') for t in ts)
        gates['planned_execution_complete']=bool(ts) and all(t['execution']['status']=='COMPLETED' and t.get('performance',{}).get('actual_steps')==t['parameters'].get('steps') for t in ts)
        gates['units']=all(t['parameters'].get('locked_units')==u.as_dict() for t in ts)
        gates['density']=bool(equilibrium) and all(t.get('mass_conserved') and abs(t['density_si']/case['density_kg_m3']-1)<=c['proposed_additional']['global_density_relative_error'] for t in equilibrium)
        material=[t for t in ts if t['kind'] in ('equilibrium','flow')]
        gates['temperature']=bool(material) and all(t.get('temperature_status')=='PASS_PROPOSED' for t in material)
        gates['stationarity']=bool(material) and all(t.get('stationarity_status')=='PASS_PROPOSED' for t in material)
        target=case['kinematic_viscosity_m2_s'];width=tol['viscosity_target_relative_error']
        def viscosity_pass(t):
            v=t.get('viscosity',{});ci=v.get('ci95_si');nu=v.get('nu_si')
            return bool(ci and nu and v.get('sampling_status')=='SUFFICIENT' and ci[0]>=(1-width)*target and ci[1]<=(1+width)*target
                        and max(ci[1]-nu,nu-ci[0])/nu<=tol['viscosity_ci_relative_halfwidth'] and v.get('relative_rms',1)<=tol['profile_relative_rms'])
        gates['viscosity']=bool(base) and all(viscosity_pass(t) for t in flows)
        gates['step_and_force_sensitivity']=sens['status']=='PASS_PROPOSED'
        gates['mechanical_pressure_definition']=bool(material) and all(t.get('pressure_channels_verified') for t in material)
        gates['EOS_statistics']=eos['status']=='MEASURED_LOCAL_RESPONSE'
        gates['EOS_temperature']=eos.get('thermal_target_status')=='PASS_PROPOSED'
        gates['pressure_range']=eos.get('target_pressure_coverage_status')=='COVERED'
        gates['density_adjustment']=bool(eos['points']) and all(abs(p['rho_star']/u.to_star(case['density_kg_m3'],'mass_density')-1)<=c['proposed_additional']['eos_density_change_max'] for p in eos['points'])
        cs=eos.get('sound_speed_star');mach=None
        if cs and flows:
            mach=max(max(abs(v) for v in t.get('viscosity',{}).get('measured_profile_star',[math.inf]))/cs for t in flows)
        gates['measured_EOS_Mach']=mach is not None and math.isfinite(mach) and mach<=tol['mach_max']
        if method=='SDPD':
            gates['local_density_fluctuation']=bool(material) and all(t.get('kernel_density',{}).get('status')=='PASS_PROPOSED' for t in material)
            gates['kernel_density_observation']=bool(material) and all(t.get('kernel_CPU_audit',{}).get('status')=='PASS' for t in material)
            gates['temperature_bin_refinement']=bool(flows) and all(t.get('temperature_refinement',{}).get('status')=='PASS_PROPOSED' for t in flows)
        # DPD kernel density is not its force-law variable; equilibrium mechanical EOS still must be tested.
        reasons=[k for k,passed in gates.items() if not passed]
        cost_only=all(t['kind']=='probe' for t in ts)
        status='COST_ONLY' if cost_only else ('QUALIFIED_PROPOSED' if not reasons else 'NOT_QUALIFIED')
        rows.append({'method':method,'candidate_id':cid,'comparison_group':group,'status':status,'gates':gates,'failed_or_unverified':reasons,
          'base_flow_task':base['task_id'] if base else None,'base_viscosity':base.get('viscosity') if base else None,
          'equilibrium_temperature_star':[t.get('temperature',{}).get('mean') for t in equilibrium],
          'flow_temperature_star':base.get('temperature',{}).get('mean') if base else None,
          'Mach_from_measured_EOS':mach if mach is not None and math.isfinite(mach) else None,'eos':eos,'sensitivity':sens,
          'time_to_qualified_solution_s':sum(t['execution']['elapsed_monotonic_s'] for t in ts) if status=='QUALIFIED_PROPOSED' else None,
          'not_validated':['bulk viscosity','wall no-slip and permeability','kernel truncation at walls','RBC membrane coupling','open pressure reservoirs']})
    return rows


def costs(tasks,assessments):
    by={t['task_id']:t for t in tasks if t['source_category']=='MEASURED'}
    names=['sdpd_probe','dpd_cost','dpd_cost_density'];matched=[by[n] for n in names if n in by]
    keys=['N','domain_star','rc_star','precision','sample_interval_physical_s','actual_steps']
    comparable=len(matched)==3 and all(all(t['performance'][k]==matched[0]['performance'][k] for k in keys) for t in matched)
    records=[{'method':t['method'],'candidate_id':t['candidate_id'],'task_id':t['task_id'],**t['performance']} for t in tasks if 'performance' in t]
    estimate=None
    if comparable:
        a=by['dpd_cost']['performance']['synchronized_steady_compute_s_per_step'];b=by['dpd_cost_density']['performance']['synchronized_steady_compute_s_per_step']
        estimate={'extra_s_per_step':b-a,'relative_to_plain_DPD':b/a-1,
                  'scope':'Matched passive Density path includes intermediate halos, channel preservation and scheduling; not an isolated CUDA-kernel timer. Single short trial, no statistical speed CI.'}
    qualified={a['method'] for a in assessments if a['status']=='QUALIFIED_PROPOSED'}
    return {'execution_cost':{'matched_short_control_status':'COMPARABLE_EXECUTION_SETTINGS' if comparable else 'INCOMPLETE_OR_UNMATCHED',
               'controls':[t['task_id'] for t in matched],'density_increment':estimate,'runs':records,
               'scope':'These parameters have different material properties. Short-cost controls do not establish target-liquid speed or large-vessel speedup.'},
            'target_liquid_efficiency':{'ranking_allowed':qualified=={'DPD','SDPD'},'ranking':None,
               'time_to_qualified_solution_s':{a['method']+'/'+a['candidate_id']:a['time_to_qualified_solution_s'] for a in assessments},
               'reason':'Only compare target-liquid efficiency after both methods reach similar physical properties and uncertainty. No speed winner is inferred from failed materials.'}}


def export_results(c):
    case,mapping,u=frozen_targets(c);old=verify_historical(c);campaign=output(Path(c['runs_root'])/c['campaign_id'])
    new=read_json(campaign/'budget_ledger.json') if (campaign/'budget_ledger.json').exists() else {'attempts':[]}
    if any(a['status']=='RUNNING' for a in new['attempts']):raise RuntimeError('EXPERIMENT_STILL_RUNNING')
    sources=[PROJECT_ROOT/'py_scripts/fluid_comparison',PROJECT_ROOT/'py_scripts/fluid_physics']
    code={str(p):sha256_file(p) for root in sources for p in root.glob('*.py')}
    hashes={a['directory']:sha256_file(Path(a['directory'])/'output_sha256.json') for a in old['attempts']+new['attempts']}
    identity={'config':c['_config_sha256'],'historical_manifest':c['historical_manifest_sha256'],'runs':hashes,'implementation':code}
    validation_ref=campaign/'validation_record.json'
    validation=read_json(validation_ref) if validation_ref.exists() else {'status':'NOT_RECORDED'}
    if validation_ref.exists():
        if sha256_file(Path(validation['log_path']))!=validation['log_sha256']:raise ValueError('CPU_TEST_LOG_HASH_MISMATCH')
        identity['CPU_validation']=validation
    key=fingerprint(identity);d=output(Path(c['data_root'])/('result_'+key[:16]))
    if d.exists():verify_package(d);return d
    oldc=yaml.safe_load((Path(c['historical_result'])/'fluid_physics.yaml').read_text())
    oldc['legacy_windows']=c['legacy_windows'];oldc['proposed_additional']=c['proposed_additional']
    tasks=[]
    for historical,ledger,config in [(True,old,oldc),(False,new,c)]:
        for a in ledger['attempts']:
            if (Path(a['directory'])/'output_sha256.json').exists():
                tasks.append(analyze_run(a['directory'],config,u,case,historical=historical))
            else:
                spec=a['cache_identity']['task'];cand=a['cache_identity']['candidate']
                tasks.append({'method':cand['method'],'candidate_id':cand['id'],'comparison_group':spec['comparison_group'],'test_kind':spec['test_kind'],
                              'task_id':spec['id'],'kind':spec['kind'],'parameters':{},'status':'FAILED_NO_COMPLETE_MANIFEST','execution':{'status':a['status'],'elapsed_monotonic_s':a['charged_s']},
                              'source_category':'MEASURED','directory':a['directory']})
    # Each failed attempt remains visible; only the last successful/planned record per identity is evidence.
    latest={ (t['source_category'],t['method'],t['candidate_id'],t['comparison_group'],t['task_id']):t for t in tasks }
    assessments=assess(list(latest.values()),case,u,c);selection=choose(assessments)
    budget=shared_budget_state(campaign,c);budget={k:v for k,v in budget.items() if k!='pool'}
    budget['historical_used_s']=sum(a['charged_s'] for a in old['attempts']);budget['new_used_s']=sum(a['charged_s'] for a in new['attempts'])
    cost=costs(list(latest.values()),assessments)
    executed={a['task_id'] for a in new['attempts']}
    summary={'schema_version':1,'status':'PARTIAL' if selection is None else 'PROPOSED','human_review':'PENDING',
             'tasks':tasks,'candidate_assessments':assessments,'selection':selection,'costs':cost,'shared_budget':budget,
             'unexecuted_planned_tasks':[t for t in c['experiments'] if t['id'] not in executed],
             'old_analysis_changes':['Candidate-local EOS and sensitivity selection','EOS bulk modulus uses actual mass density','Nominal pressure coverage separated from conservative baseline/endpoint CI coverage; uncertainty never enlarges measured range','Legacy windows read from frozen acceptance mapping','Warmup checks both flow and thermal stationarity','Timestamp tolerance uses the actual six-significant-digit serialization quantum','Historic UTC annotation corrected in new report only'],
             'scope':'Periodic homogeneous pure liquid only. No real vessel, SDF, RBC, bubbles, adhesion, walls or open-boundary simulation.'}
    contract=native_contract()
    verified=[t['task_id'] for t in tasks if t['method']=='SDPD' and t.get('pressure_channels_verified') and t.get('kernel_CPU_audit',{}).get('status')=='PASS']
    contract['runtime_verified_tasks']=verified
    if verified:contract['status']='SOURCE_AND_NATIVE_OBSERVATIONS_VERIFIED'
    d.mkdir(parents=True,exist_ok=False)
    for name,value in [('physical_case',case),('unit_mapping',mapping),('native_sdpd_contract',contract),('CPU_validation',validation),
                       ('validation_checks',{'scientific_status':summary['status'],'human_review':'PENDING','CPU_tests':validation.get('status','NOT_RECORDED'),
                         'raw_hashes_verified':True,'JSON_nonfinite_forbidden':True,'package_roundtrip':'VERIFIED_BY_EXPORT_READER','GPU_budget_within_authorization':budget['total_charged_or_reserved_s']<=budget['limit_s'],
                         'all_experiments_visible':True,'source_and_runtime_contract':contract['status'],'selection_requires_all_independent_gates':True}),
                       ('comparison_summary',summary),('selected_parameters',{'selection':selection,'status':summary['status'],'human_review':'PENDING'}),
                       ('provenance',{'identity':identity,'created_at':now(),'historical_reference':c['historical_result'],'new_campaign':str(campaign),'environment':environment_identity(),
                                      'measurement_categories':['INPUT','DERIVED','MEASURED','FITTED','HISTORICAL_REFERENCE','PROPOSED','UNVERIFIED']}),
                       ('resolved_config',c),('measurement_methods',{'viscosity':'Inherited periodic force/mass/bin-average analytic fit, independent measured slope and block CI',
                         'temperature':'New flow uses instantaneous half-width bins and unbiased degrees of freedom; equilibrium subtracts COM. Old coarse-bin residual estimated separately.',
                         'pressure':'Native total pair virial plus COM-subtracted kinetic momentum flux once, divided by V. Only equilibrium mechanical pressure enters EOS evidence.',
                         'local_EOS':'Input formula evaluated at observed m*d_i; separately tagged and never used as an independent mechanical measurement.',
                         'sampling':'Each independent initialization analysed separately; inherited warmup search 20-50%, physical relaxation minimum, ACF and nonoverlapping blocks. Insufficient windows never PASS.',
                         'performance':'Completion-guaranteed synchronized compute blocks, host/output times, setup and full charged process wall time; device memory includes other applications.'})]:
        write_json(d/(name+'.json'),value)
    shutil.copyfile(c['_config_path'],d/'fluid_model_comparison.yaml')
    with (d/'comparison.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['method','candidate','status','input_mu_star','measured_nu_m2_s','target_error','CI95_low','CI95_high','flow_temperature_star','time_to_qualified_s'])
        for a in assessments:
            v=a['base_viscosity'] or {};ci=v.get('ci95_si') or [None,None]
            w.writerow([a['method'],a['candidate_id'],a['status'],v.get('input_mu_star'),v.get('nu_si'),v.get('target_relative_error'),*ci,a['flow_temperature_star'],a['time_to_qualified_solution_s']])
    (d/'report_zh.md').write_text(chinese_report(case,mapping,summary))
    write_json(d/'package_sha256.json',{p.name:sha256_file(p) for p in d.iterdir() if p.is_file()});verify_package(d)
    return d


def chinese_report(case,mapping,s):
    def fmt(x):return '未得到有效结果' if x is None else f'{x:.7g}'
    b=s['shared_budget'];lines=['# DPD 与原生 SDPD 纯液体对照','',f"状态：**{s['status']}**；selection={json.dumps(s['selection'],ensure_ascii=False)}；用户人工验收 PENDING。",
      '',f"固定目标：ρ={case['density_kg_m3']} kg/m³，ν={case['kinematic_viscosity_m2_s']} m²/s，μ={case['dynamic_viscosity_pa_s']} Pa·s，T={case['temperature_K']} K。",
      f"原保存单位 L0={mapping['L0']:.17g} m，M0={mapping['M0']:.17g} kg，t0={mapping['t0']:.17g} s；目标 ν*={mapping['required_nu_star']:.12g}、μ*={mapping['required_mu_star']:.12g}。没有重新定义时间单位。",
      '', '## 材料判定','', '| 模型/候选 | SDPD 输入 μ* | 实测 ν (m²/s) | 95% CI | 绝对相对误差 | 流动态 T/T目标 | 结论 |','|---|---:|---:|---|---:|---:|---|']
    for a in s['candidate_assessments']:
        v=a['base_viscosity'] or {};ci=v.get('ci95_si')
        ci_text='无有效 CI' if ci is None else '['+', '.join(fmt(x) for x in ci)+']'
        error='未测' if v.get('target_relative_error') is None else f"{100*v['target_relative_error']:.3f}%"
        mu='不适用' if v.get('input_mu_star') is None else fmt(v['input_mu_star'])
        lines.append(f"|{a['method']}/{a['candidate_id']}|{mu}|{fmt(v.get('nu_si'))}|{ci_text}|{error}|{fmt(a['flow_temperature_star'])}|{a['status']}|")
    for a in s['candidate_assessments']:
        equilibrium=', '.join(fmt(x) for x in a['equilibrium_temperature_star']) or '未测'
        lines += ['',f"### {a['method']} / {a['candidate_id']}",f"未通过或缺证：{', '.join(a['failed_or_unverified']) or '无'}。",
          f"平衡温度 T*：{equilibrium}；步长/驱动：{a['sensitivity']['status']}；机械 EOS：{a['eos']['status']}；压力覆盖：{a['eos'].get('target_pressure_coverage_status','UNVERIFIED')}。"]
        v=a['base_viscosity'] or {}
        if v:
            lines.append(f"基准流动黏度有效统计块 {v.get('block_valid_count',0)} 个；统计状态 {v.get('sampling_status','UNVERIFIED')}；剖面相对 RMS {fmt(v.get('relative_rms'))}。有效块不足时不报告有效 CI，也不把均值接近目标解释为通过。")
    lines += ['', '## 实际执行与费用','',f"旧账本已用 {b['historical_used_s']:.6f} s；本轮新增 {b['new_used_s']:.6f} s；合计 {b['total_charged_or_reserved_s']:.6f} / {b['limit_s']:.0f} s；剩余 {b['remaining_s']:.6f} s。",
      '', '下表的统计段时长是所选记录覆盖的物理区间，统计充分性另判；不是保证独立、合格的采样时长。',
      '', '| 任务 | 来源 | 墙钟 s | 实际物理时长 s | 选取统计段 s | 每步总成本 s | 退出状态 |','|---|---|---:|---:|---:|---:|---|']
    for t in s['tasks']:
        p=t.get('performance',{})
        lines.append(f"|{t['task_id']}|{t['source_category']}|{fmt(t['execution']['elapsed_monotonic_s'])}|{fmt(p.get('actual_physical_duration_s'))}|{fmt(p.get('statistics_physical_duration_s'))}|{fmt(p.get('total_wall_s_per_step'))}|{t['execution']['status']}|")
    execution=s['costs']['execution_cost'];increment=execution['density_increment']
    lines += ['',f"同观测短成本对照：{execution['matched_short_control_status']}。以下为实际短任务设置的执行成本，物性不同，不是合格目标液体的效率排名。",
      '', '| 短成本任务 | 同步计算 ms/step | 完整墙钟 s/物理 µs | 初始化 s | 采样/传输/输出 s |','|---|---:|---:|---:|---:|']
    for t in s['tasks']:
        if t['source_category']!='MEASURED' or t['task_id'] not in execution['controls']:continue
        p=t['performance']
        lines.append(f"|{t['task_id']}|{fmt(1000*p['synchronized_steady_compute_s_per_step'])}|{fmt(1e-6*p['wall_s_per_physical_s'])}|{fmt(p['worker_setup_wall_s'])}|{fmt(p['sampling_transfer_output_s'])}|")
    if increment:
        lines += ['',f"DPD 增加被动 Density 路径后，同步计算成本增加 {1000*increment['extra_s_per_step']:.6f} ms/step（{100*increment['relative_to_plain_DPD']:.3f}%）。该差值包含密度通信、通道保存和调度，不是孤立 CUDA 核函数计时。"]
    lines += ['',
      '完整计时、初始化、平衡段、同步计算、采样传输输出、主机可用内存与设备显存位于 comparison_summary.json。设备显存含其他应用；短测试每步成本没有重复实验的置信区间。',
      f"目标液体效率排名允许：{'是' if s['costs']['target_liquid_efficiency']['ranking_allowed'] else '否'}。未达到合格材料时，达到合格解的耗时为 null；不能按低黏度历史基线的速度宣布血流加速。",
      '', '## 测量定义与局限',
      'SDPD 输入的是动力黏度，实测运动黏度来自独立速度剖面拟合。模拟点、目标解析线和拟合线分别保存。',
      '全盒 N/V、N*m/V、核数密度及进入 EOS 的 m*d_i 分开保存；Wendland C2 包含自贡献。少量原始快照由 CPU 独立重算密度，只作观测核查，CPU 耗时单列。',
      '机械压力采用原生配对应力的迹，加一次动能动量通量并除以体积；不再添加 EOS 压力。EOS 局部预测与机械压力的差异保留，不能当作已校准的绝对热力学压力。',
      '旧三个表压差原值保留。每个候选以本身的实测平衡机械压力作为测试零表压参考，不套用旧 LBM 数值偏置或大气压；移动零点不扩大已测压力区间。',
      '压力均值名义覆盖与包含不确定性的覆盖分开。保守核查要求端点压力的向内 CI 边界覆盖“基准压力加表压差”的完整基准 CI；各点是单独 95% CI，不宣称联合 95% 覆盖。',
      '原生 Stats/virial 时间差由该版本两个回调和每次实际 CSV 对应验证。核密度快照取 beforeIntegration，和同相位位置一起保存；其时间与步末速度相差一个实际 dt。',
      '局部材料统计门槛与旧血管短/长窗口分别报告。本轮不是空间收敛测试；相同间距不保证相同空间误差。没有将独立初始化拼成连续轨迹。',
      '', '## 旧数据重新审查','旧温度偏高同时发生于平衡态和流动态，因此不能只归咎于未扣流速。各流动记录另给基于实测剖面的有限箱宽剪切残余估计；没有通过减去任意速度场强行通过温度检查。',
      '旧高耗散候选的时间步/驱动力与 EOS 证据独立判定，不借用低耗散候选。旧 CSV、配置、报告、UTC 字段和账本全部保留；分析修正只写入本结果。',
      '', '## 后续能力边界','未运行真实血管、SDF、RBC、微泡、黏附、管壁或开放压力出口。后续需要验证近壁核截断、冻结壁粒子、膜内外液体/膜节点耦合，以及双向压力储液区和有限截面压力测量。普通 DPD/SDPD 的体积黏度未独立匹配。',
      '', '## 下一步','依据本轮实际失败项先解决初始化/温度和统计充分性，再验证黏度及机械 EOS；保留当前未合格结论，不自动选失败候选为生产默认。']
    return '\n'.join(lines)+'\n'
